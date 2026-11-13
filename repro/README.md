# 旋转连续性复现与边界扫描

所有程序均为纯 C，依赖库头文件与 `build/libcglm.so`。

## 构建

```sh
mkdir -p build && cd build
cmake -DCGLM_USE_TEST=ON .. && make -j4
cd ..
cc -O2 -std=c99 -I include repro/repro.c          -L build -lcglm -lm -o build/repro
cc -O2 -std=c99 -I include repro/boundary_scan.c  -L build -lcglm -lm -o build/boundary_scan
cc -O2 -std=c99 -I include repro/pole_check.c     -L build -lcglm -lm -o build/pole_check
cc -O2 -std=c99 -I include repro/long_arc.c       -L build -lcglm -lm -o build/long_arc
```

## 程序说明

- `repro.c` — 主复现用例，4 组场景：
  1. 单轴连续 5 圈、逐帧归一化后的关键帧 slerp（含跨 ±π 缝与近对极回退）；
  2. 多圈欧拉角 `euler -> matrix -> euler` 往返；
  3. 俯仰扫过 ±π/2 万向节锁点，校验重建姿态连续且精确锁点不丢滚转；
  4. 729 点网格 `euler -> quat -> matrix -> euler` 往返误差。
- `boundary_scan.c` — 边界扫描，输出 `boundary_scan.csv`（归档于 `artifacts/`），
  含 ±π 缝、±π/2 极点、精确锁点、4913 点多圈网格和 slerp 跨缝统计。
- `pole_check.c` — 用双精度参考实现对比，证明近极点误差按 1/距离 线性增长，
  属于欧拉角固有条件数而非跳变。
- `long_arc.c` — 短弧/长弧 slerp 中点与端点不变量。

## 运行

```sh
LD_LIBRARY_PATH=build ./build/repro
LD_LIBRARY_PATH=build ./build/boundary_scan   # 生成 boundary_scan.csv
LD_LIBRARY_PATH=build ./build/pole_check
LD_LIBRARY_PATH=build ./build/long_arc
```
