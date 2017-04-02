#!/bin/bash

#./maplesat  -ccmin-mode=0 -never-gc -always-restart -lsr-num -lsr-out=lsr $1 | tee genlsr
#./maplesat  -ccmin-mode=0 -never-gc -always-restart -lsr-cert-cls-out=cert -lsr-in=lsr $1 | tee gencert
#./maplesat  -ccmin-mode=0 -never-gc -always-restart -lsr-cert-cls-in=cert -lsr-in=lsr $1 | tee checkcert

./maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -always-restart -conf-side-lsr -verb=0 $1
./maplesat -lsr-in=lsr  -lsr-cert-cls-out=clauses -lsr-num -ccmin-mode=2 -always-restart -conf-side-lsr -verb=0 $1
./maplesat -lsr-in=lsr  -lsr-cert-cls-in=clauses -never-gc -lsr-num -ccmin-mode=2 -always-restart -conf-side-lsr -verb=0 $1
