################################################################################################
function GS_Periodic!(du, u ,p::Array{Float64, 1},t::Float64) # Works with square and rect grids
  f, k, D₁, D₂, dx, dy, N = p

  let N = Int(N)
    @floop for j in 2:N-1, i in 2:N-1
      du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
    end

    @floop for j in 2:N-1, i in 2:N-1
      du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
    end

    let i = 1
      @floop for j in 2:N-1
        du[i,j,1] = D₁*(1/dx^2*(u[N,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                    -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
        du[i,j,2] = D₂*(1/dx^2*(u[N,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                    u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
        end
    end

    let i = N
      @floop for j in 2:N-1
        du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                    -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
        du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                    u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
      end
    end

    let j = 1
      @floop for i in 2:N-1
        du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,N,1] - 2u[i,j,1])) +
                    -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
        du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,N,2] - 2u[i,j,2])) +
                    u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
      end
    end

    let j = N
      @floop for i in 2:N-1
        du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                    -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
        du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                    u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
      end
    end

    begin
      i = 1; j = 1
      du[i,j,1] = D₁*(1/dx^2*(u[N,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,N,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
      du[i,j,2] = D₂*(1/dx^2*(u[N,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,N,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]

      i = 1; j = N
      du[i,j,1] = D₁*(1/dx^2*(u[N,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
      du[i,j,2] = D₂*(1/dx^2*(u[N,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]

      i = N; j = 1
      du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,N,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
      du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,N,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]

      i = N; j = N
      du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
      du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
    end
  end

  nothing
end

function GS_Neumann0!(du, u, p::Array{Float64, 1},t::Float64) # Works with square and rect grids
  f, k, D₁, D₂, dx, dy, N = p

  let  N = Int(N)
    @floop for j in 2:N-1, i in 2:N-1
      du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
    end

    @floop for j in 2:N-1, i in 2:N-1
      du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
    end

    let i = 1
      @floop for j in 2:N-1
        du[i,j,1] = D₁*(1/dx^2*(2u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                    -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
        du[i,j,2] = D₂*(1/dx^2*(2u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                    u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
      end
    end

    let i = N
      @floop for j in 2:N-1
        du[i,j,1] = D₁*(1/dx^2*(2u[i-1,j,1] - 2u[i,j,1])+ 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                    -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
        du[i,j,2] = D₂*(1/dx^2*(2u[i-1,j,2] - 2u[i,j,2])+ 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                    u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
      end
    end

    let j = 1
      @floop for i in 2:N-1
        du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j+1,1] - 2u[i,j,1])) +
                    -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
        du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j+1,2] - 2u[i,j,2])) +
                    u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
      end
    end

    let j = N
      @floop for i in 2:N-1
        du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j-1,1] - 2u[i,j,1])) +
                    -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
        du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j-1,2] - 2u[i,j,2])) +
                    u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
      end
    end

    begin
      i = 1; j = 1
      du[i,j,1] = D₁*(1/dx^2*(2u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j+1,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
      du[i,j,2] = D₂*(1/dx^2*(2u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j+1,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]

      i = 1; j = N
      du[i,j,1] = D₁*(1/dx^2*(2u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j-1,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
      du[i,j,2] = D₂*(1/dx^2*(2u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j-1,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]

      i = N; j = 1
      du[i,j,1] = D₁*(1/dx^2*(2u[i-1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j+1,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
      du[i,j,2] = D₂*(1/dx^2*(2u[i-1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j+1,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]

      i = N; j = N
      du[i,j,1] = D₁*(1/dx^2*(2u[i-1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j-1,1] - 2u[i,j,1])) +
                  -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
      du[i,j,2] = D₂*(1/dx^2*(2u[i-1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j-1,2] - 2u[i,j,2])) +
                  u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
    end
  end
  nothing 
end