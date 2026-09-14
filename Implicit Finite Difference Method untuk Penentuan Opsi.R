# Implicit Finite Difference Method untuk Penentuan Opsi
# Berdasarkan Black-Scholes-Merton Model

# Parameter
S0 <- 145    # Harga saham saat ini
K <- 140     # Harga exercise
r <- 0.0025  # Tingkat bunga bebas risiko
T <- 1       # Masa kadaluwarsa
sigma <- 0.387443624 # Volatilitas
M <- 100   # Jumlah step
N <- 100   # Jumlah step dari harga saham

# Ukuran waktu dan harga step
dt <- T/M
ds <- S0/N

# Buat kisi harga dan waktu saham
s <- seq(0, N*ds, by=ds)
t <- seq(0, M*dt, by=dt)

# Menginisialisasi matriks harga opsi
V <- matrix(0, nrow=N+1, ncol=M+1)

# Tetapkan kondisi batas
V[,M+1] <- pmax(s-K,0)  # Saat jatuh tempo
V[1,] <- 0              # Batas bawah harga saham
V[N+1,] <- N*ds-K*exp(-r*(T-t))  # Batas atas harga saham

# Hitung koefisien dari matriks tridiagonal
a <- 0.5*dt*(sigma^2*s^2-r*s)
b <- 1+dt*(sigma^2*s^2+r)
c <- -0.5*dt*(sigma^2*s^2+r*s)

# Fungsi untuk mencari solusi dari sistem matriks tridiagonal dengan Algoritma Thomas  

solve.TDMA <- function(a, b, c, d) {
  n <- length(d)
  cprime <- rep(0, n-1)
  dprime <- rep(0, n)
  x <- rep(0, n)
  cprime[1] <- c[1]/b[1]
  dprime[1] <- d[1]/b[1]
  for (i in 2:(n-1)) {
    cprime[i] <- c[i]/(b[i]-a[i]*cprime[i-1])
  }
  for (i in 2:n) {
    dprime[i] <- (d[i]-a[i]*dprime[i-1])/(b[i]-a[i]*cprime[i-1])
  }
  x[n] <- dprime[n]
  for (i in (n-1):1) {
    x[i] <- dprime[i]-cprime[i]*x[i+1]
  }
  return(x)
}

# Solusi dari sistem tridiagonal dengan algoritma Thomas 

for (j in M:1) {
  d <- V[,j+1]
  d[1] <- d[1] - a[1]*V[1,j]
  d[N+1] <- d[N+1] - c[N+1]*V[N+1,j]
  V[,j] <- solve.TDMA(a, b, c, d)
}

# Interpolasi harga opsi pada harga saham awal/inisial

price <- approxfun(s, V[,1], method="linear")(S0)

# Munculkan harga opsi

cat("Option price:", price)

