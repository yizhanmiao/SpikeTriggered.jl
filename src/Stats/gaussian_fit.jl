# function Gaussian2DSimplex(X, Y)
#     (param) -> begin
#         (α, σ_a, σ_b, θ, x0, y0, bias) = length(param) == 6 ? [param; 0] : param

#         a = cos(θ^2)/(2 * σ_a^2) + sin(θ^2)/(2 * σ_b^2)
#         b = sin(2*θ)/(4 * σ_a^2) - sin(2*θ)/(4 * σ_b^2)
#         c = sin(θ^2)/(2 * σ_a^2) + cos(θ^2)/(2 * σ_b^2)
#         _deltaX = X[:, 1] .- x0
#         _deltaY = X[:, 2] .- y0
#         _raw = @. α * exp(-(a * _deltaX^2 + b * _deltaX * _deltaY + c * _deltaY^2)) + bias

#         return mean(abs2, _raw.-Y)
#     end
# end

# function _make_ellipse(a, b, θ, x0, y0; σ=1, step=100)
#     t = range(0; stop=2*pi, length=step)
#     ellipse_x_r = a .* cos.(t)
# 	ellipse_y_r = b .* sin.(t)
#     R = [cos(-θ) sin(-θ); -sin(-θ) cos(-θ)] .* σ
#     r_ellipse = [ellipse_x_r ellipse_y_r] * R
#     x = x0 .+ r_ellipse[:,1]
# 	y = y0 .+ r_ellipse[:,2]
#     (x, y)
# end

# #NOTE: visualization with Makie
# # function ellipse!(scene, args...; σ=1, step=100, kwargs...)
# #     lines!(scene, _make_ellipse(args...; σ, step)...; kwargs...)
# # end

# function Gaussian2DSimplexInit(x::AbstractMatrix, y::AbstractVector{T}; bias=true) where {T}
#     param = T[
#         maximum(y);
#         1;
#         1;
#         0.1;
#         mean(x, weights(y), 1)[:];
#     ]
#     bias && push!(param, mean(y))
#     return param
# end

# #TODO: docs
# function srf_gaussian_fit(image::AbstractMatrix; bias=true, optim_options=(;), full=false)
#     X, Y = image2dataset(image)
#     loss = Gaussian2DSimplex(X, Y)
#     param_init = Gaussian2DSimplexInit(X, Y; bias)
#     rez = optimize(loss, param_init, NelderMead(), OptimOptions(;optim_options...))
#     return (full ? rez : rez.minimizer)
# end

# #TODO: migrate `srf` submodule here
