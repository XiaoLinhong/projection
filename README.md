# projection

支持三种安装方式。

使用MAKE安装
```
make run
```

使用CMAKE管理
```
mkdir build; cd build

FC=ifort cmake ../ # 正常编译
FC=ifort cmake  -DCMAKE_BUILD_TYPE=Debug ../ # 调试

make
```

使用fortran包管理器(fpm)
```
fpm build
```
