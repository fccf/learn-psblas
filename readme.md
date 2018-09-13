
PSBLAS(parallel sparse BLAS) 主要用于并行迭代求解大型稀疏线性代数方程组。 该程序用Fortran2003按照面向对象思想编写。

<!-- more -->
## 安装
### 下载
``` bash
git clone https://github.com/sfilippone/psblas3.git
```
### 安装
按照说明进行安装，基本上就是以下命令
```
./configure
make
make install
make clean
```
要预先安装blas，lapack，metis等库。

有时configure找不到blas等库文件，需要手动指定
```
./configure FCOPT=-O3 CCOPT=-O3 FDEFINES='-DHAVE_BLAS'
```

有时 make 需要 sudo 权限。

### 测试
使用FoBis编译测试文件，定义变量 **$PSB_DIR**，然后修改以下的变量就可以使用psblas的所有库了。
```
$PSB_DIR     = xxx


mpi             = True
include         = $PSB_DIR/modules $PSB_DIR/include
libs            = -L$PSB_DIR/lib
ext_libs        =  blas lapack metis psb_util psb_krylov psb_prec psb_base
```
## 概念
### 点的分类
  - 内部点(inner points)
    只依赖于同一区域的点；

  - 边界点(boundary points)
    依赖于其它区域的点；

  - 光晕点(halo points)(ghost points)
    一个区域的光晕点是该区域的边界点所依赖的其他区域的点；一个区域的边界点通常是其它区域的光晕点。

对于一个区域形成的矩阵，行数为 $I+B$ 行，列数为 $I+B+H$ 列，其中 $I$ 为内部点数目， $B$为边界点数目，$H$ 为光晕点数目。

### 索引空间(index_space)
建立索引空间是计算的第一步，将索引分配到不同的进程，建立全局与局部的对应关系。
+ 全局索引(global indices)
  - 全局行
  - 全局列
+ 局部索引(local indices)
  - 局部行
  - 局部光晕点 -> 局部列
  - 局部全局对照表

不同的离散方法有不同的索引空间，要确定索引空间必须确定网格的稀疏模式。

### 内置的稀疏矩阵存储格式
```
- coo          坐标存储(ia(nnz),ja(nnz),a(nnz))
- csr          压缩行存储(ia(nr+1),ja(nnz),a(nnz))
- csc          压缩列存储
```

### 矩阵的状态
```
- build        分配空间之后，组装之前
- assembled    组装后，之后可以进行矩阵计算。
- update       对相同稀疏模式矩阵值进行更新
```

## 命名原则
- 所有的subroutine 和data types 以psb_ 开头
- 所有自定义的类以 \_type 结尾
- 自定义的数据精度
  + s 单精度实数
  + d 双精度实数
  + c 单精度复数
  + z 双精度复数
- 所有常数用 _ 结尾
- 所有用户API 命名方式为 psb_xxname
  + xx主要有
    - ge: 密集数据
    - sp: 稀疏数据
    - cd: 数据交换描述符
  + name主要有
    - init -> initialize
    - all  -> allocate
    - ins  -> insert
    - asb  -> assemble
    - prec -> precondition
    - bld  -> build

## 主要API
### 空间分配API

- 通信描述符分配psb_cdall

```
call psb_cdall(ictxt,desc_a,info,mg=mg,parts=parts)
call psb_cdall(ictxt,desc_a,info,vg=vg,[mg=mg,flag=flag])
call psb_cdall(ictxt,desc_a,info,vl=vl,[nl=nl,globalcheck=.true.,lidx=lidx])
call psb_cdall(ictxt,desc_a,info,nl=nl)
call psb_cdall(ictxt,desc_a,info,mg=mg,repl=.true.)

integer(psb_ipk_)     :: ictxt   ! communication context
type(psb_desc_type)   :: desc_a  ! communication descriptor
integer(psb_ipk_)     :: info    ! Error code
integer(psb_ipk_)     :: mg      ! the (global) number of rows of the problem
external              :: parts   ! the partitioning subroutine
integer(psb_ipk_)     :: vg(mg)  ! each index(i) is allocated to process vg(i)
integer(psb_ipk_)     :: nl      ! number of local indices
integer(psb_ipk_)     :: vl(nl)  ! global indices of current process
logical               :: repl    ! build a replicated index space

```
cd (communication descriptor)存储进程之间的通信信息。通信描述符分配将离散的点分配到不同的核，第一次调用后进入build state，此时可以在不同的进程中添加通信。通过psb_cdins或psb_spins添加边，即定义整个离散的模式或稀疏矩阵的稀疏模式。然后通过psb_cdasb或psb_spasb计算通信信息，即定义halo points等，进入assembled state。此时就可以进行并行的矩阵计算，内部通信会自动进行。

- 稀疏矩阵分配psb_spall
```
call psb_spall(a, desc_a, info, nnz)
type(psb_dspmat_type) :: a      ! sparse matrix
type(psb_desc_type)   :: desc_a ! communication descriptor
integer(psb_ipk_)     :: info   ! Error code
integer(psb_ipk_)     :: nnz    ! An estimate of the number of nonzeroes of local

```


- 密集矩阵或向量分配psb_geall

```
call psb_geall(x, desc_a, info, n, lb)
type(psb_d_vect_type) :: x      ! dense vector or matrix
type(psb_desc_type)   :: desc_a ! communication descriptor
integer(psb_ipk_)     :: info   ! Error code
integer(psb_ipk_)     :: n      ! The number of columns (matrix)
integer(psb_ipk_)     :: lb     ! The lower bound for the column index (matrix)
```

### 赋值API
- 通信描述符赋值
```
call psb_cdins(nz, ia, ja, desc_a, info [,ila,jla])
call psb_cdins(nz,ja,desc,info[,jla,mask,lidx])

```
- 稀疏矩阵赋值
```
call psb_spins(nz, ia, ja, val, a, desc_a, info [,local])
integer(psb_ipk_)     :: nz
integer(psb_ipk_)     :: ia(nz)
integer(psb_ipk_)     :: ja(nz)
integer(psb_dpk_)     :: val(nz)
type(psb_dspmat_type) :: a      ! sparse matrix
type(psb_desc_type)   :: desc_a ! communication descriptor
integer(psb_ipk_)     :: info   ! Error code
```

- 密集矩阵或向量赋值
```
call psb_geins(m, irw, val, x, desc_a, info [,dupl,local])
type(psb_d_vect_type) :: x      ! dense vector or matrix
type(psb_desc_type)   :: desc_a ! communication descriptor
integer(psb_ipk_)     :: info   ! Error code
```
赋值API定义整个稀疏矩阵的稀疏模式以及稀疏矩阵的值，

### 组装API
- 通信描述符组装
```
call psb_cdasb(desc_a, info)
```
- 稀疏矩阵组装
```
call psb_spasb(a, desc_a, info, afmt, upd, dupl, mold)
type(psb_dspmat_type) :: a      ! sparse matrix
type(psb_desc_type)   :: desc_a ! communication descriptor
integer(psb_ipk_)     :: info   ! Error code

```

- 密集矩阵或向量组装
```
call psb_geasb(x, desc_a, info, mold)
```

组装完成之后，整个cd定义完成，即定义好给个区域的halo points，定义好了每个进程的发送和接收的信息，进入assembled状态，之后就可以进行并行的矩阵计算以及方程求解，不用担心内部的通信问题。


### 并行环境API
```
- psb_init     并行环境初始化
- psb_info     返回并行环境信息
- psb_exit     退出并行环境
- psb_abort    中断计算
- psb_barrier  暂时中断
- psb_wtime    返回时间
```

### 数据管理API

### 计算API

## 计算流程
```
- 并行环境初始化 psb_init
- 索引空间初始化 psb_cdall
- 稀疏矩阵初始化 psb_spall
- 矩阵向量初始化 psb_geall
- 各进程系数计算 psb_spins, psb_geins
- 矩阵组装      psb_cdasb, psb_spasb, psb_geasb
- 预处理        psb_precset, psb_precbld
- 迭代求解      psb_krylov
```
