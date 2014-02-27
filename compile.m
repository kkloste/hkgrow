
if ismac
% mex -O -largeArrayDims hkseed_mex.cpp
mex -O -largeArrayDims hkgrow_mex.cpp
% mex -O -largeArrayDims hkgrow_degs_mex.cpp
else
% mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims hkseed_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims hkgrow_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims hkgrow_degs_mex.cpp
end