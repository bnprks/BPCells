#// Copyright 2025 BPCells contributors
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#include <condition_variable>
#include <mutex>
#include <thread>

#include "DenseTransposeMultiply.h"

// Overall approach:
//   - Load from input to in-memory CSC sparse matrix chunks
//   - Multiply each chunk with its transpose and update the dense output
//       - Doing an explicit transpose ensures that we have better cache locality
//         when updating the output, as we can apply all outputs to a single column at a time
//   - Use multi-threading to split up the multiply work, and buffer chunks so that
//     data loading can overlap with multiply computations

// Internal helper classes and functions
namespace {

/**
 * Represent a compressed sparse matrix in memory (CSC or CSR format)
 * - val, idx (length >= # non-zeros)
 * - ptr (length # cols/rows + 1)
 * Unusual features in this code:
 * - ptr might temporarily be modified during other algorithms
 * - val and idx might have extra items beyond the number of non-zeros, which should be left
 *   untouched by other algorithms
 */
struct CSparseMat {
    std::vector<double> val;
    std::vector<uint32_t> idx;
    std::vector<size_t> ptr;
};

/**
 * Initialize the out.ptr vector
 * Preconditions:
 *  - `in` is a valid CSparseMat, though in.idx and in.val might have extra values at the end which
 * should be left in place.
 *  - out.ptr is sized correctly
 * Postconditions:
 *  - out.ptr is a valid pointer index for the transposed compressed storage view
 */
inline void init_transpose_ptr(const CSparseMat &in, CSparseMat &out) {
    for (auto &x : out.ptr) {
        x = 0;
    }
    for (size_t i = 0; i < in.ptr.back(); i++) {
        out.ptr[in.idx[i] + 1]++;
    }
    for (size_t i = 1; i < out.ptr.size(); i++) {
        out.ptr[i] += out.ptr[i - 1];
    }
}

/**
 * Transpose compressed sparse matrix storage.
 * Preconditions:
 *   - `in` is a valid CSparseMat, though in.idx and in.val might have extra values at the end which
 * should be left in place. (same for out.idx and out.val having extra usage beyond the number of
 * nonzeros)
 *   - `out.ptr` has been initialized by `init_transpose_ptr`
 * Postconditions:
 *   - `out` holds a transposed copy of `in`
 */
inline void transpose(const CSparseMat &in, CSparseMat &out) {
    if (out.idx.size() < in.ptr.back()) out.idx.resize(in.ptr.back());
    if (out.val.size() < in.ptr.back()) out.val.resize(in.ptr.back());

    for (size_t c = 0; c < in.ptr.size() - 1; c++) {
        for (size_t i = in.ptr[c]; i < in.ptr[c + 1]; i++) {
            size_t i_out = out.ptr[in.idx[i]];
            out.val[i_out] = in.val[i];
            out.idx[i_out] = c;
            out.ptr[in.idx[i]]++;
        }
    }

    // Shift row_ptr values to give the correct pointers again
    std::memmove(
        out.ptr.data() + 1,
        out.ptr.data(),
        (out.ptr.size() - 1) * sizeof(decltype(out.ptr)::value_type)
    );
    out.ptr[0] = 0;
}

/**
 * Perform a subset of the sparse multiply operations
 * @param trans Row-major sparse storage of the column chunk
 * @param m Col-major sparse storage of the column chunk
 * @param ret Output dense results matrix
 * @param task_id ID of the task (< n_tasks)
 * @param n_tasks How many total tasks the multiply is being split into
 *
 * To try to balance work among tasks, each task calculates output entries of
 * every n_tasks columns.
 */
void chunk_multiply_worker(
    const CSparseMat &trans,
    const CSparseMat &m,
    Eigen::MatrixXd &ret,
    size_t task_id,
    size_t n_tasks
) {
    for (size_t r = task_id; r < trans.ptr.size() - 1; r += n_tasks) {
        double *out = &ret(0, r);
        for (size_t i = trans.ptr[r]; i < trans.ptr[r + 1]; i++) {
            uint32_t c = trans.idx[i];
            double vi = trans.val[i];
            for (size_t j = m.ptr[c]; j < m.ptr[c + 1]; j++) {
                if (m.idx[j] >= r) {
                    out[m.idx[j]] += vi * m.val[j];
                }
            }
        }
    }
}

/**
 * Load CSparseMat chunks from a MatrixLoader
 */
class SparseChunkLoader {
  private:
    std::unique_ptr<BPCells::MatrixLoader<double>> loader_;
    const size_t buf_items_;
    bool finished_matrix_ = false;
    bool finished_column_ = true;

    CSparseMat mat_;

  public:
    SparseChunkLoader(std::unique_ptr<BPCells::MatrixLoader<double>> &&mat, size_t buf_items)
        : loader_(std::move(mat))
        , buf_items_(buf_items) {

        mat_.ptr.push_back(0);
        mat_.val.reserve(buf_items);
        mat_.idx.reserve(buf_items);
    }

    // Return a reference to the CSparseMat chunk
    CSparseMat &mat();

    /** Load a new compressed sparse chunk of data
     *  @return true if new data was loaded, false if the matrix loader is finished
     */
    bool next();

    // Swap the contents of the internal matrix with another CSparseMat,
    // preserving any leftovers at the end of the internal matrix and overwriting `other`
    void swap_mat(CSparseMat &other);
};

// Allow a leader thread to coordinate with worker threads
class WorkerBarrier {
  private:
    const size_t workers_;
    size_t waiting_ = 0;
    size_t epoch_ = 0;
    bool finished = false;
    std::mutex m_;
    std::condition_variable worker_cv_, leader_cv_;

  public:
    WorkerBarrier(size_t workers) : workers_(workers) {}

    // Call from leader thread. Resumes execution for all worker threads without blocking
    void start_workers() {
        {
            std::unique_lock<std::mutex> lk(m_);
            leader_cv_.wait(lk, [this] { return waiting_ == workers_; });
            epoch_++;
            waiting_ = 0;
        }
        worker_cv_.notify_all();
    }

    // Call from leader thread. Blocks until all workers have reached a wait point
    void wait_for_workers() {
        std::unique_lock<std::mutex> lk(m_);
        leader_cv_.wait(lk, [this] { return waiting_ == workers_; });
    }

    // Call from a worker thread after performing work.
    // Will pause execution until leader thread next calls start_workers() or finish().
    // Returns false once no more work is left to be done (i.e. finish has been called)
    bool wait_for_leader() {
        std::unique_lock<std::mutex> lk(m_);
        size_t run = epoch_;
        waiting_++;
        if (waiting_ == workers_) {
            lk.unlock();
            leader_cv_.notify_one();
            lk.lock();
        }
        worker_cv_.wait(lk, [this, run] { return epoch_ != run || finished; });
        return !finished;
    }

    // Call from leader thread. Signals to worker threads to finish at their next wait point.
    void finish() {
        {
            std::unique_lock<std::mutex> lk(m_);
            finished = true;
        }
        worker_cv_.notify_all();
    }

    const size_t workers() { return workers_; }
};

} // namespace

Eigen::MatrixXd dense_transpose_multiply(
    std::unique_ptr<BPCells::MatrixLoader<double>> &&mat,
    size_t buffer_bytes,
    size_t threads,
    std::atomic<bool> *user_interrupt
) {
    if (buffer_bytes < 48 * mat->rows()) {
        throw std::runtime_error(
            "dense_transpose_multiply: buffer_bytes must be at least 24 * mat.rows()"
        );
    }

    threads = std::max(threads, (size_t) 1);

    std::vector<std::thread> thread_vec;
    thread_vec.reserve(threads - 1);

    Eigen::MatrixXd ret(mat->rows(), mat->rows());
    ret.setZero();

    size_t buf_items = buffer_bytes / (12 * 4);
    if (threads == 1) buf_items *= 2;

    CSparseMat trans, trans_next, m;
    trans.ptr.resize(mat->rows() + 1);
    trans_next.ptr.resize(mat->rows() + 1);

    SparseChunkLoader l(std::move(mat), buf_items);
    CSparseMat &m_next = l.mat();

    // Handle single-threaded case standalone
    if (threads <= 1) {
        while (l.next()) {
            init_transpose_ptr(m_next, trans_next);
            transpose(m_next, trans_next);
            transpose(trans_next, m_next);
            chunk_multiply_worker(trans_next, m_next, ret, 0, 1);
            if (user_interrupt != NULL && *user_interrupt) break;
        }

        ret = ret.selfadjointView<Eigen::Lower>();
        return ret;
    }

    // Set up multiplication worker threads:
    // pull from a shared "queue" with more tasks than threads in case of uneven task lengths
    WorkerBarrier b(threads - 1);
    const size_t tasks = std::min(threads * 10, trans.ptr.size() - 1);
    std::atomic<size_t> task_id(0);
    auto run_worker = [&]() {
        while (b.wait_for_leader()) {
            while (true) {
                size_t id = task_id.fetch_add(1);
                if (id >= tasks) break;
                chunk_multiply_worker(trans, m, ret, id, tasks);
            }
        }
    };
    for (size_t i = 0; i < threads - 1; i++) {
        thread_vec.emplace_back(run_worker);
    }

    // Load first batch of data
    bool finished = !l.next();
    init_transpose_ptr(m_next, trans_next);
    transpose(m_next, trans_next);
    // Second transpose ensures are column-major chunk has nonzeros orderd by row, which is common
    // but not guaranteed straight from the MatrixLoader interface.
    transpose(trans_next, m_next); 

    while (!finished) {
        l.swap_mat(m);
        std::swap(trans_next, trans);

        // Run workers while loading next data
        task_id = 0;
        b.start_workers();

        finished = !l.next();
        init_transpose_ptr(m_next, trans_next);
        transpose(m_next, trans_next);
        transpose(trans_next, m_next); // See above comment for why second transpose is needed

        // If the main thread has finished a load before the worker threads
        // are done, then do some chunk_multiply work tasks
        while (true) {
            size_t id = task_id.fetch_add(1);
            if (id >= tasks) break;
            chunk_multiply_worker(trans, m, ret, id, tasks);
        }

        if (user_interrupt != NULL && *user_interrupt) break;
        b.wait_for_workers();
    }
    b.wait_for_workers();
    l.swap_mat(m);
    std::swap(trans_next, trans);
    b.start_workers();
    b.wait_for_workers();
    b.finish();

    for (auto &t : thread_vec) {
        t.join();
    }

    return ret.selfadjointView<Eigen::Lower>();
}

namespace {

CSparseMat &SparseChunkLoader::mat() { return mat_; }

bool SparseChunkLoader::next() {
    if (finished_matrix_) return false;

    // Move any partial column data to the front of the buffers
    size_t leftover = mat_.val.size() - mat_.ptr.back();
    std::memmove(
        mat_.val.data(),
        mat_.val.data() + mat_.ptr.back(),
        leftover * sizeof(decltype(mat_.val)::value_type)
    );
    mat_.val.resize(leftover);
    std::memmove(
        mat_.idx.data(),
        mat_.idx.data() + mat_.ptr.back(),
        leftover * sizeof(decltype(mat_.idx)::value_type)
    );
    mat_.idx.resize(leftover);

    mat_.ptr.clear();
    mat_.ptr.push_back(0);

    uint32_t first_col = 0;
    while (true) {
        if (finished_column_ && !loader_->nextCol()) {
            finished_matrix_ = true;
            return true;
        }

        while (!finished_column_ || loader_->load()) {
            finished_column_ = true;
            uint32_t cap = loader_->capacity();
            if (cap + mat_.val.size() > buf_items_) {
                finished_column_ = false;
                return true;
            }
            mat_.val.insert(mat_.val.end(), loader_->valData(), loader_->valData() + cap);
            mat_.idx.insert(mat_.idx.end(), loader_->rowData(), loader_->rowData() + cap);
        }
        finished_column_ = true;
        if (mat_.ptr.size() == 1) first_col = loader_->currentCol();
        // In case loader_->nextCol() ever skipped a column (I don't think the current interfaces do
        // that, but being safe)
        if (loader_->currentCol() + 1 > mat_.ptr.size() + first_col) {
            mat_.ptr.insert(
                mat_.ptr.end(), (loader_->currentCol() - first_col) - (mat_.ptr.size() - 1), 0
            );
        }
        mat_.ptr.push_back(mat_.val.size());
    }
}

void SparseChunkLoader::swap_mat(CSparseMat &other) {
    // Transfer any partial column data from mat_ to other as needed.
    size_t leftover = mat_.val.size() - mat_.ptr.back();
    other.val.resize(leftover);
    other.idx.resize(leftover);
    other.ptr.resize(1);
    other.ptr[0] = 0;
    std::memmove(
        other.val.data(),
        mat_.val.data() + mat_.ptr.back(),
        leftover * sizeof(decltype(mat_.val)::value_type)
    );
    std::memmove(
        other.idx.data(),
        mat_.idx.data() + mat_.ptr.back(),
        leftover * sizeof(decltype(mat_.idx)::value_type)
    );
    mat_.val.resize(mat_.ptr.back());
    mat_.idx.resize(mat_.ptr.back());

    // Perform the swap
    std::swap(mat_, other);
}

} // namespace
