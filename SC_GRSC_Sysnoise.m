
function  [Sys_noise , iter]     =     SC_GRSC_Sysnoise( Img_list,Opts)

b              =   Opts.win;
[h, w,img_num]     =   size(Img_list);
N              =   h-b+1;
M              =   w-b+1;

nsig           =   Opts.nSig;
m              =   Opts.nblk;
gamma          =   Opts.gamma;

Out_Put        =   Img_list;
Sys_noise=zeros(h, w);
X0=zeros(b*b,N*M,img_num);
for nn=1:img_num
    X0(:,:,nn) =Im2Patch( Img_list(:,:,nn), Opts ); 
end

for iter = 1 : Opts.Iter
    disp(iter)
    dif= Sys_noise;
    
    vd=nsig^2-(mean(mean(dif.^2)));
    if iter==1                    
        Opts.nSig =nsig; 
    else
        Opts.nSig =sqrt(abs(vd))*Opts.lamada;          
    end
    Ys_NEW        =     zeros( b*b,N*M );  
    W             =     zeros( b*b,N*M  );
    
    E_Img  	=  zeros(h,w);
    W_Img 	=  zeros(h,w);
    for nn=1:img_num
        
        disp(nn)
        Out_Put(:,:,nn)=Img_list(:,:,nn)-Sys_noise + gamma*Sys_noise;
        X =Im2Patch( Out_Put(:,:,nn), Opts );  
        NL_mat =Block_matching( Out_Put(:,:,nn), Opts);
        K             =     size(NL_mat,2);

        for  i  =  1 : K  

            % Get Nonlocal Similar patches from noisy image...
           A             =      X(:, NL_mat(:, i));
           A_Nim=X0(:, NL_mat(:, i),nn);


           TMP_NEW            =      GRSC( double(A), Opts.c1, Opts.nSig, Opts.eps,Opts.hp,double(A_Nim));
           Ys_NEW(:, NL_mat(1:m,i))    =   Ys_NEW(:, NL_mat(1:m,i)) + TMP_NEW;
           W(:, NL_mat(1:m,i))     =   W(:, NL_mat(1:m,i)) + 1;
    
        end
        [E_Img,W_Img]=  Patch2Im(  E_Img,W_Img,Ys_NEW, W, b, h, w );
    end
    
    Sys_noise  =  E_Img./(W_Img+eps);
    
    if  iter >1
      
       dif=norm(abs(Sys_noise) - abs(dif),'fro')/norm(abs(dif), 'fro');
       
       if dif<Opts.errr_or
           
           break;
       end
       
   end
        
end
