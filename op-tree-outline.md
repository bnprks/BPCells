

Rough idea for how the interfaces might get used
```cpp
// Define an op with its type-checking data
auto op = OpBuilder("ShiftCoords")
    .set_inputs({OpType::Fragment})
    .set_output(OpType::Fragment)
    .add_param(ParamBuilder<int32_t>("shift_start").scalar().default(0))
    .add_param<int32_t>("shift_end", Scalar(0))
    .set_executor(&ShiftCoords::execute);

registry.addOp(op);

// Execute an op
registry.execute<MatrixLoader>(OpTree);

// Get analysis data (a.k.a computed metadata) from an operation
// (Uncertain API)
registry.fragmentInfo(OpTree)
    .name()
    .chromosomeType();
```

OpBuilder:
    - set_inputs()