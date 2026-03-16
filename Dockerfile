FROM julia:1.11
LABEL maintainer="Morgan Askins <maskins@berkeley.edu>"

RUN julia -e 'using Pkg; \
              Pkg.add(url="https://github.com/MorganAskins/Batman.jl.git"); \
              using Batman'
