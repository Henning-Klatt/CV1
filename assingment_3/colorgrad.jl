# Author: Christian Jarvers
#
# Copyright 2020 Christian Jarvers (MIT License)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and # to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of 
# the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
# THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.

module ColorGrad

using Images

export colorgrad

"""
    colorgrad(angles, magnitudes, threshold)

Generate a color image where hue is given by `angles` and intensity is given by `magnitudes`.
`magnitudes` smaller than `threshold` are suppressed.
"""
function colorgrad(angles, magnitudes, threshold)
    # threshold and normalize magnitudes
    m = max.(magnitudes, threshold)
    m = (m .- minimum(m)) / (maximum(m) - minimum(m))
    # generate HSV image
    img = colorview(HSV, rad2deg.(angles), ones(size(angles)), m)
    return img
end
    
end # module ColorGrad
