
if ismac
mex -O -largeArrayDims hkgrow_mex.cpp
else
mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims hkgrow_mex.cpp
end