!
        integer, parameter :: ns_itrmax = 10	! max no. of iter in trace estimation
        integer, parameter :: nst0      =  5	! initial guess for no. of samples

        real    trace_tol                       ! relative tolerance
        logical do_rndtrace                     ! when .t. will calc trace via rnd trace algorithm
        logical trace_corr                      ! when .t. calculate trace of correlation matrix
        common / rndtr01 / do_rndtrace, trace_tol, trace_corr
