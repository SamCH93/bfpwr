## Source locally so tinytest's masked options() restores the caller's settings
## when the test file exits, including early exits and errors. This profile is
## only for behavior tests; default and numerical-accuracy tests stay unchanged.
options(bfpwr.ngrid = 1000, bfpwr.tail.nquad = 128)
