mex -largeArrayDims -O triangleclustersgreedy_mex.cc
mex -largeArrayDims -O triangleclusters_mex.cc
if ismac
    mex -largeArrayDims greedyclustergrow.cc
    mex -largeArrayDims pprgrow_mex.cc
    mex -largeArrayDims cutcond_mex.cc
else
    mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims greedyclustergrow.cc
    mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims pprgrow_mex.cc
    mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims cutcond_mex.cc
end