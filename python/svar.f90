module model_t

  use, intrinsic :: iso_fortran_env, only: wp => real64
  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
  use fortress_prior_t, only: fortress_abstract_prior, M_PI, lognorpdf, logigpdf
  use fortress_random_t, only: fortress_random
  use fortress_linalg, only: cholesky, Kronecker, determinant
  use fortress_util, only: read_array_from_file
  use logbeta, only : betaln, gamln

  implicit none

  type, public, extends(fortress_abstract_prior) :: minnesota_prior

     integer  :: p = {p}, constant = {cons:d}, ny = {ny}
     integer  :: nA = {nA}, nF = {nF}
     
     real(wp) :: hyper_phistar({nF}), hyper_Omega_inv({nF}/{ny},{nF}/{ny}), hyper_iw_Psi({ny},{ny})
     integer :: hyper_iw_nu 

   contains

     procedure rvs
     procedure logpdf
     procedure para_to_sigma_phi
     procedure para_to_A_F


  end type minnesota_prior

  interface minnesota_prior
     module procedure new_prior
  end interface minnesota_prior



  type, public, extends(fortress_abstract_bayesian_model) :: model

     integer :: p = {p}, constant = {cons:d}
     integer :: nA = {nA}, nF = {nF}
     integer :: likT 

     real(wp), allocatable :: YY_lik(:,:), XX_lik(:,:), proxy(:)

   contains
     procedure :: lik

  end type model

  interface model
     module procedure new_model
  end interface model

contains

  type(model) function new_model() result(self)

    character(len=144) :: datafile, name
    
    integer :: i,j

    !self%npara = self%nA + self%nF
    datafile = 'data.txt'
    name = '{name}'
    call self%construct_abstract_bayesian_model(name, datafile, self%nA+self%nF+{nextra_para}, {ny}, {T})
    self%likT = self%T - self%p

    allocate(self%prior,source=minnesota_prior())

    allocate(self%YY_lik(self%likT, self%nobs), &
         self%XX_lik(self%likT, self%nobs*self%p+self%constant), self%proxy(self%likT))

    self%XX_lik = 1.0_wp
    do i = 1, self%likT
       self%YY_lik(i,:) = self%YY(:,self%p+i)
       do j = 1, self%p
          self%XX_lik(i,(j-1)*self%nobs+1:j*self%nobs) = self%YY(:,i-j+self%P)
       end do
    end do

    self%T = self%likT

    call read_array_from_file('proxy.txt', self%proxy)
  end function new_model

  real(wp) function lik(self, para, T) result(l)

    class(model), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    integer, intent(in), optional :: T

    integer :: use_T, ny, nF, ipiv(self%nobs), p, constant, TT, info, i

    real(wp) :: A(self%nobs,self%nobs), F(self%nF/self%nobs, self%nobs), work(3*self%nobs)

    real(wp) :: likvec(self%likT), A0(self%nobs, self%nobs), AXiAXip(self%nobs,self%nobs)
    real(wp) :: resid(self%likT, self%nobs), beta, signu
    
    use_T = self%likT
    ny = self%nobs
    nF = self%nF
    p = self%p
    TT= self%likT
    constant = self%constant
    if (present(T)) use_T = T

    associate(prior => self%prior)
      select type(prior)
      class is (minnesota_prior)
         call prior%para_to_A_F(para, A, F)
      class default
         print*,'prior misspecified'
         stop
      end select
    end associate
       

    ! ! zero out financial shock
    ! A(2,3:5) = 0.0_wp
    ! A(3,4:5) = 0.0_WP
    ! A(4,5) = 0.0_wp

    call dgemm('n','n',TT,ny,ny,1.0_wp,self%YY_lik,TT,A,ny,0.0_wp,resid,TT)
    call dgemm('n','n',TT,ny,ny*p+constant,-1.0_wp,self%XX_lik,TT,F,ny*p+constant,1.0_wp,resid,TT)

    beta = para(self%nA+self%nF+1)
    signu = {signu}



    resid(:,1) = ( self%proxy - beta*resid(:,1) ) / signu
    l = -TT/2.0_wp*log(2.0_wp*M_PI*signu**2) - 0.5_wp*sum(resid(:,1)**2)

    
  end function lik

  type(minnesota_prior) function new_prior() result(pr)

    pr%npara = pr%nA + pr%nF + {nextra_para}
    
    call read_array_from_file('phistar.txt', pr%hyper_phistar)
    call read_array_from_file('Omega_inv.txt',pr%hyper_Omega_inv)
    call read_array_from_file('iw_Psi.txt',pr%hyper_iw_Psi)
    
    pr%hyper_iw_nu = {nu}

  end function new_prior


  function rvs(self, nsim, seed, rng) result(parasim)

    class(minnesota_prior), intent(inout) :: self

    integer, intent(in) :: nsim
    integer, optional :: seed
    type(fortress_random), optional, intent(inout) :: rng

    real(wp) :: parasim(self%npara, nsim)

    ! pass -- read in from file
  end function rvs
    

  real(wp) function logpdf(self, para) result(lpdf)

    class(minnesota_prior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    real(wp) :: sigma(self%ny,self%ny), phi(self%nF/self%ny, self%ny), phi_vec(self%nF)

    real(wp) :: mvn_covar(self%nF, self%nF)

    integer :: ny, nF, info, p, i

    ny = self%ny
    nF = self%nF
    p = self%p
    lpdf = 0.0_wp


    call self%para_to_sigma_phi(para, sigma, phi)

    lpdf = logiwishpdf(self%hyper_iw_nu, self%hyper_iw_Psi, sigma, self%ny)

    call Kronecker(ny,ny,nF/ny,nF/ny,nF,nF,0.0_wp,1.0_wp,sigma,self%hyper_Omega_inv,mvn_covar)
    call cholesky(mvn_covar, info)

    do i = 1, self%ny
       phi_vec((i-1)*self%nF/self%ny+1:i*self%nF/self%ny) = phi(:,i)
    end do
    lpdf = lpdf + mvnormal_pdf(phi_vec, self%hyper_phistar, mvn_covar)

    ! ARRW
    lpdf = lpdf + (ny+1.0_wp)/2.0_wp*log(2.0_wp) + (2*ny+ny*p+2)/2.0_wp*log(determinant_gen(sigma, ny))

    lpdf = lpdf + lognorpdf(para(self%nA+self%nF+1), 0.0_wp, 0.1_wp)

    if (self%npara==self%nA+self%nF+2) then
       if (para(self%nA+self%nF+2) < 0.0_wp) then
          lpdf = -100000000000.0_wp
          return
       end if

       lpdf = lpdf + logigpdf(para(self%nA+self%nF+2), 0.02_wp, 2.0_wp)
    end if

  end function logpdf
    

  subroutine para_to_A_F(self, para, A, F)

    class(minnesota_prior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    real(wp), intent(out) :: A(self%ny, self%ny), F(self%nF/self%ny, self%ny)

    {assign_para}

  end subroutine para_to_A_F

  subroutine para_to_sigma_phi(self, para, sigma, phi)

    class(minnesota_prior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    real(wp) :: A(self%ny, self%ny), F(self%nF/self%ny, self%ny)
    real(wp), intent(out) :: sigma(self%ny,self%ny), phi(self%nF/self%ny, self%ny)

    real(wp) :: A0i(self%ny, self%ny), A0(self%ny, self%ny), work(3*self%ny)
    integer :: k, ind0, ny, nF, info, ipiv(self%ny), col

    call self%para_to_A_F(para, A, F)

    ny = self%ny
    nF = self%nF
    A0 = A
    A0i = A0
    call dgetrf(ny,ny,A0i,ny,ipiv,info)
    call dgetri(ny,A0i,ny,ipiv,work,3*ny,info)
    col = nF/ny
    phi = 0.0_wp
    call dgemm('n','n',col,ny,ny,1.0_wp,F,col,A0i,ny,0.0_wp,phi,col)
  
    call dgemm('n','t',ny,ny,ny,1.0_wp,A0,ny,A0,ny,0.0_wp,sigma,ny)
    call dgetrf(ny,ny,sigma,ny,ipiv,info)
    call dgetri(ny,sigma,ny,ipiv,work,3*ny,info)
 
 

  end subroutine para_to_sigma_phi


  function logiwishpdf(nu,S,X,n)

    integer, intent(in) :: nu, n
    real(wp), intent(in) :: S(n,n), X(n,n)
    real(wp) :: logiwishpdf
    integer :: i, ipiv(n), info
    real(wp) :: kap, iX(n,n), tr, work(n), ret(n,n),x0

    iX = X
    call dgetrf(n, n, iX, n, ipiv, info)
    call dgetri(n, iX, n, ipiv, work, n, info)
    call dgemm('n','n',n,n,n,1.0_wp,iX,n,S,n,0.0_wp,ret,n)

    tr = 0.0_wp
    kap = 0.0_wp
    x0 = (nu*1.0_wp)/2.0_wp

    do i = 1,n
       tr = tr + ret(i,i)
       kap = kap + gamln(x0+(1-i)*0.5_wp)!gamln(0.5_wp*(nu*1.0_wp-i*1.0_wp-1.0_wp))
    end do

    kap = -kap - 0.5_wp*nu*n*log(2.0_wp) - 0.25_wp*n*(n-1)*log(M_PI)


    logiwishpdf = kap + 0.5_wp*nu*log(determinant(S,n)) &
         - 0.5_wp*(nu+n+1)*log(determinant(X,n)) -0.5_wp*tr

  end function logiwishpdf


  function mvnormal_pdf(x, mu, chol_sigma) result(logq)
    !! Computes the log of the n-dimensional multivariate normal pdf at x
    !! with mean mu and variance = chol_sigma*chol_sigma'.
    !!
    real(wp), intent(in) :: x(:), mu(:), chol_sigma(:,:)
    real(wp), external :: ddot

    real(wp) :: logq, lds
    real(wp) :: a(size(x, 1)), det_sigma

    integer :: n, i
    n = size(x, 1)

    if (n /= size(mu, 1)) then
       print*, 'mvnormal pdf, size error'
       stop
    endif

    det_sigma = 1.0_wp
    lds = 0.0_wp
    do i = 1, n
       lds = lds + log(chol_sigma(i,i))
       det_sigma = det_sigma*chol_sigma(i,i)

    end do
    det_sigma = det_sigma**2
    a = x - mu
    call dtrsv('l','n', 'n', n, chol_sigma, n, a, 1)

    logq = -n*0.5_wp*log(2.0_wp*3.14159_wp) - lds - 0.5*ddot(n, a, 1, a, 1)

  end function mvnormal_pdf

function determinant_gen(mat,n) result(det)
 
  integer, intent(in) :: n
  real(wp), intent(in) :: mat(n,n)
 
  integer :: i, info, ipiv(n)
 
  real(wp) :: det
  real(wp) :: sgn, matcp(n,n)
 
  det = 1.0_wp
  sgn = 1.0_wp
  matcp = mat
  call dgetrf(n,n,matcp,n,ipiv,info)
 
  do i = 1,n
     det = det*matcp(i,i)
     if (ipiv(i)/=i) then
        sgn = -sgn
     end if
  end do
 
  det = sgn*det
 
end function determinant_gen
 

end module model_t
