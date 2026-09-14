# Input data
S0 <- 145    # Harga saham saat ini
K <- 140     # Harga exercise
r <- 0.0025  # Tingkat bunga bebas risiko
T <- 1       # Masa kadaluwarsa
sigma <- 0.387443624 # Volatilitas

# Grid parameters
M <- 100    # Jumlah grid points
N <- 10000  # Jumlah time steps
dt <- T/N   # Increment waktu
dx <- log(S0/K)/M   # Increment harga saham

# Inisialisasi matriks untuk hasil
v <- matrix(0, nrow = M+1, ncol = N+1)

# Inisialisasi kondisi batas
v[,1] <- pmax(S0*exp(seq(0,-M*dx,-dx))-K,0)
v[1,] <- 0
v[M+1,] <- S0*exp(-r*T)-K*exp(-r*(T-seq(0,N)*dt))

# Inisialisasi matriks tridiagonal
alpha <- 0.25*dt*(sigma^2*(1:M)^2-r*(1:M))
beta <- -dt*0.5*(sigma^2*(1:M)^2+r)
gamma <- 0.25*dt*(sigma^2*(1:M)^2+r*(1:M))

A <- diag(-beta,1) + diag(alpha[-1],1) + diag(gamma[-(M-1)],1)
B <- diag(beta,1) - diag(alpha[-1],1) - diag(gamma[-(M-1)],1)


# Metode Crank-Nicolson
for (i in 2:(N+1)) {
  # Solusi sistem linier Ax = Bv_i-1
  x <- solve(A, B %*% v[-c(1,M+1),i-1])
  
  # Update harga opsi
  v[-c(1,M+1),i] <- x
}

# Harga opsi
v0 <- v[(M+1)/2,N+1]

# Print hasil
cat("Harga opsi Eropa:", round(v0,2))
