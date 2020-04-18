pyHilbertSort
==============================================================================================
Bin Wang (binwang.0213@gmail.com), Craft and Hawkins Department of Petroleum Engineering, Louisiana State University, Baton Rouge, US

A light-weight hilbert curve sorting library for point and mesh, core algorithm is modified from Zoltan @ http://www.cs.sandia.gov/Zoltan/dev_html/dev_hsfc.html

<p align="center">
  <img src = "https://github.com/BinWang0213/HilbertSort/blob/master/examples/OrigOrder.png" height="250">
  <img src = "https://github.com/BinWang0213/HilbertSort/blob/master/examples/Reorder.png" height="250">
</p>

## Prerequisites

**On Unix (Linux, OS X)**

* A compiler with C++11 support
* CMake >= 2.8.12
* Eigen 3.3.7
* pybind11

**On Windows**

* Visual Studio 2015 (required for all Python versions, see notes below)
* CMake >= 3.1
* Eigen 3.3.7
* pybind11


## Installation

```bash
pip install pyHilbertSort
```

Or build from source: Just clone this repository and pip install. Note the `--recursive` option which is
needed for the pybind11 submodule:

```bash
conda install -c anaconda cmake
conda install -c conda-forge eigen

git clone --recursive https://github.com/BinWang0213/HilbertSort.git
cd ./HilbertSort
python setup.py install 
```

The `python setup.py install` command will
invoke CMake and build the code as specified in `CMakeLists.txt`.

## Usage

```python
import pyHilbertSort as m
pts_sort,idx_sort=m.hilbertSort(3, np.random.rand(5,3))
print(pts_sort,idx_sort)
```

## Special notes for Windows

**Compiler requirements**

Pybind11 requires a C++11 compliant compiler, i.e Visual Studio 2015 on Windows.
This applies to all Python versions, including 2.7. Unlike regular C extension
modules, it's perfectly fine to compile a pybind11 module with a VS version newer
than the target Python's VS version. See the [FAQ] for more details.

**Runtime requirements**

The Visual C++ 2015 redistributable packages are a runtime requirement for this
project. It can be found [here][vs2015_runtime]. If you use the Anaconda Python
distribution, you can add `vs2015_runtime` as a platform-dependent runtime
requirement for you package: see the `conda.recipe/meta.yaml` file in this example.


[FAQ]: http://pybind11.rtfd.io/en/latest/faq.html#working-with-ancient-visual-studio-2009-builds-on-windows
[vs2015_runtime]: https://www.microsoft.com/en-us/download/details.aspx?id=48145
