import os
env = Environment(ENV = os.environ)

env.Append(LIBS=["png","fftw3","m"])
env.Append(CXXFLAGS=["-std=c++14","-g","-Wall","-O3","-I/usr/include/libpng16"])
env.Append(LINKFLAGS=["-L/usr/local/lib"])

env.Program("stable_fluids",
            [
                "fluid.cpp", "simulate.cpp"
            ])

