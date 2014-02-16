
if ismac
mex -g -largeArrayDims hkgrow_mex.cpp
%mex -g -largeArrayDims hkpr_mex.cpp
else
mex -g -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims hkgrow_mex.cpp
%mex -g -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims hkpr_mex.cpp
end