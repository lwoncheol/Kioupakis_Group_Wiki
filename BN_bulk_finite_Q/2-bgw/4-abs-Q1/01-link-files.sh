#!/bin/bash

ln -sf ../../1-mf/4-wfn_co_full/wfn.cplx WFN_co
ln -sf ../../1-mf/4-wfnq_co_Q1/wfn.cplx WFNq_co
ln -sf ../../1-mf/5-wfn-fi-64-Q0-2/wfn.cplx WFN_fi
ln -sf ../../1-mf/6-wfnq-fi-64-Q1-2/wfn.cplx WFNq_fi
ln -sf ../3-ker-Q1/bsemat.h5 .
ln -sf ../1-eps/q0/eps0mat.h5 .
ln -sf ../1-eps/epsmat.h5 .
ln -sf ../2-sig-full/eqp1.dat eqp_co.dat
ln -sf ../2-sig-Q1/eqp1.dat eqp_co_q.dat
