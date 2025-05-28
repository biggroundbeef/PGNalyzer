from numpy import convolve, ones

def smooth(data, icut):
    return convolve(data, ones(icut) / icut, mode = "same")


