################################################################################################
function GS_Periodic!(du::Array{T,3}, u::Array{T,3} ,p::Vector{Float64},t::Float64, ex) where T <: Real
  f, k, D₁, D₂, dx, dy, M = p

  let N = Int(M)
    @floop ex for j in 1:N, i in 1:N
      @inbounds begin
        if i == 1 && j != N && j != 1 # left boundary
          du[i,j,1] = D₁*(1/dx^2*(u[N,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[N,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == N && j != N && j != 1 # right boundary
          du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif j == 1 && i != N && i != 1 # bottom boundary
          du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,N,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,N,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif j == N && i != N && i != 1 # top boundary
          du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == 1 && j == 1 # left bottom corner
          du[i,j,1] = D₁*(1/dx^2*(u[N,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,N,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[N,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,N,2] - 2u[i,j,2])) +
                    u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == 1 && j == N # left top corner
          du[i,j,1] = D₁*(1/dx^2*(u[N,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[N,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == N && j == 1 # right bottom corner
          du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,N,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,N,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == N && j == N # right top corner
          du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        else # bulk
          du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        end
      end
    end
  end
  nothing
end
##################################################################################
function GS_Neumann0!(du::Array{T,3}, u::Array{T,3} ,p::Vector{Float64},t::Float64,ex) where T <: Real
  f, k, D₁, D₂, dx, dy, M = p

  let  N = Int(M)
    @floop ex for j in 2:N-1, i in 2:N-1
        if i == 1 && j != N && j != 1 # left boundary
          du[i,j,1] = D₁*(1/dx^2*(2u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(2u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == N && j != N && j != 1 # right boundary
          du[i,j,1] = D₁*(1/dx^2*(2u[i-1,j,1] - 2u[i,j,1])+ 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(2u[i-1,j,2] - 2u[i,j,2])+ 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif j == 1 && i != N && i != 1 # bottom boundary
          du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j+1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j+1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif j == N && i != N && i != 1 # top boundary
          du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == 1 && j == 1 # left bottom corner
          du[i,j,1] = D₁*(1/dx^2*(2u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j+1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(2u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j+1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == 1 && j == N # left top corner
          du[i,j,1] = D₁*(1/dx^2*(2u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(2u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == N && j == 1 # right bottom corner
          du[i,j,1] = D₁*(1/dx^2*(2u[i-1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j+1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(2u[i-1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j+1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        elseif i == N && j == N # right top corner
          du[i,j,1] = D₁*(1/dx^2*(2u[i-1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(2u[i-1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        else # bulk
          du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                      -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
          du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                      u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        end
      end
    end
  nothing 
end
##################################################################################