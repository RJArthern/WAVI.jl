
@with_kw struct Workspace{T <: Real, N <: Integer}
nxh :: N
nyh :: N
nu :: N
nv :: N
nxu :: N = nxh + 1 
nyu :: N = nyh
nxv :: N = nxh
nyv :: N = nyh + 1
nxc :: N = nxh - 1
nyc :: N = nyh - 1
nxnyh :: N = nxh*nyh
nxnyu :: N = nxu*nyu
nxnyv :: N = nxv*nyv
nxnyc :: N = nxc*nyc
usamp :: Vector{T} = zeros(nu);                                 @assert length(usamp) == nu
vsamp :: Vector{T} = zeros(nv);                                 @assert length(vsamp) == nv
uspread :: Vector{T} = zeros(nxnyu);                            @assert length(uspread) == nxnyu
vspread :: Vector{T} = zeros(nxnyv);                            @assert length(vspread) == nxnyv
dudx :: Vector{T} = zeros(nxnyh);                               @assert length(dudx) == nxnyh
dvdy :: Vector{T} = zeros(nxnyh);                               @assert length(dvdy) == nxnyh
r_xx_strain_rate_sum :: Vector{T} = zeros(nxnyh);               @assert length(r_xx_strain_rate_sum) == nxnyh
r_yy_strain_rate_sum :: Vector{T} = zeros(nxnyh);               @assert length(r_yy_strain_rate_sum) == nxnyh
r_xx :: Vector{T} = zeros(nxnyh);                               @assert length(r_xx) == nxnyh
r_yy :: Vector{T} = zeros(nxnyh);                               @assert length(r_yy) == nxnyh
dudy_c :: Vector{T} = zeros(nxnyc);                             @assert length(dudy_c) == nxnyc
dvdx_c :: Vector{T} = zeros(nxnyc);                             @assert length(dvdx_c) == nxnyc
r_xy_strain_rate_sum_c :: Vector{T} = zeros(nxnyc);             @assert length(r_xy_strain_rate_sum_c) == nxnyc
r_xy_strain_rate_sum_crop_c :: Vector{T} = zeros(nxnyc);        @assert length(r_xy_strain_rate_sum_crop_c) == nxnyc
r_xy_strain_rate_sum :: Vector{T} = zeros(nxnyh);               @assert length(r_xy_strain_rate_sum) == nxnyh
r_xy :: Vector{T} = zeros(nxnyh);                               @assert length(r_xy) == nxnyh
r_xy_c :: Vector{T} = zeros(nxnyc);                             @assert length(r_xy_c) == nxnyc
r_xy_crop_c :: Vector{T} = zeros(nxnyc);                        @assert length(r_xy_crop_c) == nxnyc
d_rxx_dx :: Vector{T} = zeros(nxnyu);                           @assert length(d_rxx_dx) == nxnyu
d_rxy_dy :: Vector{T} = zeros(nxnyu);                           @assert length(d_rxy_dy) == nxnyu
d_ryy_dy :: Vector{T} = zeros(nxnyv);                           @assert length(d_ryy_dy) == nxnyv
d_rxy_dx :: Vector{T} = zeros(nxnyv);                           @assert length(d_rxy_dx) == nxnyv
taubx :: Vector{T} = zeros(nxnyu);                              @assert length(taubx) == nxnyu
tauby :: Vector{T} = zeros(nxnyv);                              @assert length(tauby) == nxnyv
qx :: Vector{T} = zeros(nxnyu);                                 @assert length(qx) == nxnyu
qx_crop :: Vector{T} = zeros(nxnyu);                            @assert length(qx_crop) == nxnyu
dqxdx :: Vector{T} = zeros(nxnyh);                              @assert length(dqxdx) == nxnyh
qy :: Vector{T} = zeros(nxnyv);                                 @assert length(qy) == nxnyv
qy_crop :: Vector{T} = zeros(nxnyv);                            @assert length(qy_crop) == nxnyv
dqydy :: Vector{T} = zeros(nxnyh);                              @assert length(dqydy) == nxnyh
divq :: Vector{T} = zeros(nxnyh);                               @assert length(divq) == nxnyh
extra :: Vector{T} = zeros(nxnyh);                              @assert length(extra) == nxnyh
d_extra_dx :: Vector{T} = zeros(nxnyu);                         @assert length(d_extra_dx) == nxnyu
d_extra_dy :: Vector{T} = zeros(nxnyv);                         @assert length(d_extra_dy) == nxnyv
h_d_extra_dx :: Vector{T} = zeros(nxnyu);                       @assert length(h_d_extra_dx) == nxnyu
h_d_extra_dy :: Vector{T} = zeros(nxnyv);                       @assert length(h_d_extra_dy) == nxnyv
fx :: Vector{T} = zeros(nxnyu);                                 @assert length(fx) == nxnyu
fy :: Vector{T} = zeros(nxnyv);                                 @assert length(fy) == nxnyv
fx_samp :: Vector{T} = zeros(nu);                               @assert length(fx_samp) == nu
fy_samp :: Vector{T} = zeros(nv);                               @assert length(fy_samp) == nv
opvecprod :: Vector{T} = zeros(nu+nv);                          @assert length(opvecprod) == nu + nv
end
