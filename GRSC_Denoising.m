
function  [Out_Put , iter]     =     GRSC_Denoising(noise_img, Opts)



Nim            =   noise_img;

b              =   Opts.win;

[h, w, ch]     =   size(Nim);

N              =   h-b+1;

M              =   w-b+1;
r              =   [1:N];

c              =   [1:M]; 

Out_Put        =   Nim;

gamma          =   Opts.gamma;

nsig           =   Opts.nSig;

m              =   Opts.nblk;

cnt            =   1;




for iter = 1 : Opts.Iter    
    disp(iter)
    Out_Put               =    Out_Put + gamma*(Nim - Out_Put); 
    dif=    Out_Put-Nim;   
    vd=nsig^2-(mean(mean(dif.^2)));
    [blk_arr] =    Block_matching( Out_Put, Opts);  
    
    if iter==1
        Opts.nSig =    sqrt(abs(vd)); 
    else
        Opts.nSig =    sqrt(abs(vd))*Opts.lamada;          
    end 
        
     
    X             =     Im2Patch( Out_Put, Opts );            
    Ys        =     zeros( size(X) );  
    W             =     zeros( size(X) );          
    K             =     size(blk_arr,2);
                            
    for  i  =  1 : K  

       A             =      X(:, blk_arr(:, i));         
              

       TMP  =      GRSC( double(A), Opts.c1, Opts.nSig, Opts.eps,Opts.hp);

       Ys(:, blk_arr(1:m,i))    =   Ys(:, blk_arr(1:m,i)) - TMP;
       W(:, blk_arr(1:m,i))     =   W(:, blk_arr(1:m,i)) + 1;
    
    end

    Out_Put        =   zeros(h,w);   
    im_wei        =  zeros(h,w);     
 
    k            =   0;
     for i   =  1:b
         for j  = 1:b
                k    =  k+1;
                Out_Put(r-1+i,c-1+j)  =  Out_Put(r-1+i,c-1+j) + reshape( Ys(k,:)', [N M]);
                im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + reshape( W(k,:)', [N M]);
          end
     end
     
     Out_Put            =        Out_Put./(im_wei+eps);

    
    cnt   =  cnt + 1;
        
  if  iter >1
      
       dif      =  norm(abs(Out_Put) - abs(Nim+dif),'fro')/norm(abs(Nim+dif), 'fro');
       
       if dif<Opts.errr_or
           
           break;
       end
       
  end
end
end





