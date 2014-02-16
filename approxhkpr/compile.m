
if ismac
mex -g -largeArrayDims approxhkpr_mex.cpp
else
mex -g -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims hkpr_mex.cpp
mex -g -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims approxhkpr_mex.cpp
end