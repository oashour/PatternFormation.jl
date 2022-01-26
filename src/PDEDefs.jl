################################################################################################
function GS_Periodic!(du,u,p,t) # Works with square and rect grids
  f, k, D₁, D₂, dx, dy, N = p
  N = Int(N)

  @floop for j in 2:N-1, i in 2:N-1
    du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end

  @floop for j in 2:N-1, i in 2:N-1
    du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
  end

  @floop for j in 2:N-1
    local i = 1
    du[i,j,1] = D₁*(1/dx^2*(u[N,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end
  @floop for j in 2:N-1
    local i = 1
    du[i,j,2] = D₂*(1/dx^2*(u[N,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
  end
  @floop for j in 2:N-1
    local i = N
    du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end
  @floop for j in 2:N-1
    local i = N
    du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
  end

  @floop for i in 2:N-1
    local j = 1
    du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,j+1,1] + u[i,N,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end
  @floop for i in 2:N-1
    local j = 1
    du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,j+1,2] + u[i,N,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
  end
  @floop for i in 2:N-1
    local j = N
      du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1]) + 1/dy^2*(u[i,1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end
  @floop for i in 2:N-1
    local j = N
      du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2]) + 1/dy^2*(u[i,1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
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
   nothing 
end

function GS_Neumann0!(du,u,p,t) # Works only with square grids.
  local f, k, D₁, D₂, dx, dy, M = p
  local N = Int(M)

  @floop for j in 2:N-1, i in 2:N-1
    du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end

  @floop for j in 2:N-1, i in 2:N-1
    du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
  end

  @floop for j in 2:N-1
    local i = 1
    du[i,j,1] = D₁*(1/dx^2*(2u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end
  @floop for j in 2:N-1
    local i = 1
    du[i,j,2] = D₂*(1/dx^2*(2u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
  end
  @floop for j in 2:N-1
    local i = N
    du[i,j,1] = D₁*(1/dx^2*(2u[i-1,j,1] - 2u[i,j,1])+ 1/dy^2*(u[i,j+1,1] + u[i,j-1,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end
  @floop for j in 2:N-1
    local i = N
    du[i,j,2] = D₂*(1/dx^2*(2u[i-1,j,2] - 2u[i,j,2])+ 1/dy^2*(u[i,j+1,2] + u[i,j-1,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
  end

  @floop for i in 2:N-1
    local j = 1
    du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j+1,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end
  @floop for i in 2:N-1
    local j = 1
    du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j+1,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
  end
  @floop for i in 2:N-1
    local j = N
    du[i,j,1] = D₁*(1/dx^2*(u[i-1,j,1] + u[i+1,j,1] - 2u[i,j,1])+ 1/dy^2*(2u[i,j-1,1] - 2u[i,j,1])) +
                -u[i,j,1]*u[i,j,2]^2 + f*(1-u[i,j,1])
  end
  @floop for i in 2:N-1
    local j = N
    du[i,j,2] = D₂*(1/dx^2*(u[i-1,j,2] + u[i+1,j,2] - 2u[i,j,2])+ 1/dy^2*(2u[i,j-1,2] - 2u[i,j,2])) +
                u[i,j,1]*u[i,j,2]^2 - (f+k)*u[i,j,2]
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
   nothing 
end
