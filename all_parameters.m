
function [m, CW_min, sigma, DATA, SIFS, DIFS, PIFS, ACK, delta, TUX, TUM, TCX, TCM, TXM, N, T, TsUM, TsM, CW_max]=all_parameters()
m = 5; 
CW_min = 32;
sigma=20;
DATA = 1000;  
SIFS=10;        
DIFS=50;         
PIFS=30;         
ACK=304;         
delta=2; 

TUX = DATA+SIFS+ACK+DIFS+2*delta;
TUM = DATA+SIFS+ACK+PIFS+DATA+DIFS+3*delta;
TCX = DATA+DIFS+delta;
TCM = DATA+PIFS+DATA+DIFS+2*delta;
TXM = DATA+DIFS+delta;

N=10;
T= 1e9;  

TsUM = DATA + SIFS + ACK + PIFS + 2*delta;
TsM = TXM;

CW_max = CW_min*(2^(m-1));
%==========================================================================
end