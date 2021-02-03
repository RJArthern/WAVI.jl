FROM julia:1.6

COPY . /workspaces/WAVI.jl/

RUN julia -e 'using Pkg;Pkg.add(url="/workspaces/WAVI.jl");using WAVI'

CMD ["julia"]