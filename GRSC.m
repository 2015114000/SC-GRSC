

function  Y  =   GRSC( A, B,c1, nSig,ee,A_Nim)

        
%         [U_i,SigmaY,V]=svd(full(A),'econ');    
%         PatNum=size(A,2);
%         TempC= sqrt(2)*sqrt(PatNum)*2*nSig^2;
%         %ClosedWNNM
%         temp=(SigmaY-eps).^2-4*(TempC-eps*SigmaY);
%         ind=find (temp>0);
%         svp=length(ind);
%         SigmaX=max(SigmaY(ind)-eps+sqrt(temp(ind)),0)/2; 
%         
%         SigmaX_SigmaY=zeros(size(SigmaY));
%         SigmaX_SigmaY(ind)=SigmaX./SigmaY(ind);
%         
%         
%         
% %         X=U_i(:,1:svp)*diag(SigmaX)*V(:,1:svp)' + m;   
        if nargin == 5
            A_Nim=zeros(size(A));
        end


        U_i                =              getsvd(A); % generate PCA basis
        
        A0                 =              U_i'*A;
        
        [~,m]              =              size (A0);
        
        B0                 =              U_i'*B;
            
        s0                 =              A0   -  B0;

        s0                 =              mean (s0.^2,2);

        s0                 =              max  (0, s0-nSig^2);
        
        lam                =              repmat(((c1*sqrt(2)*nSig^2)./(sqrt(s0) + ee)),[1,m]);
 
        Alpha              =              soft (A0-B0, lam)+ B0;
        

        X                  =              U_i*Alpha;
        
        Y=A_Nim-X;
        



return;