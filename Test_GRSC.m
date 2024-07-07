clear all
clc


train_path = 'traindata\';
test_path = 'testdata\';


File_train = dir(strcat(train_path,'*.tif')); 
Length_Names_train = length(File_train);
disp(Length_Names_train)

File_test = dir(strcat(test_path,'*.tif')); 
Length_Names_test = length(File_test);
disp(Length_Names_test)


%% 


image_size_h=512;
image_size_w=512;
image_size_c=1;
%% 


Sigma = 65;%system noise level
sys_noise=Sigma*randn(image_size_h,image_size_w,image_size_c);

estimated_noise=zeros(size(sys_noise));
Sys_noise_est_sum=zeros(size(sys_noise));
for k = 1 :Length_Names_train 
    disp(k)
    image_name = File_train(k).name;
    img_orig=double(imread(strcat(train_path,image_name)));

    img_noise=img_orig(:,:,1:image_size_c)+sys_noise+5*randn(size(sys_noise));

    for c=1:image_size_c
        img_noise_ch=img_noise(:,:,c);
        Opts=Opts_Set (Sigma);
        [denoised_img , iter]  =   GRSC_Denoising(img_noise_ch, Opts);


        estimated_noise(:,:,c)=img_noise_ch-(denoised_img+eps);
    end
    Sys_noise_est_sum=Sys_noise_est_sum+estimated_noise;
    
end
Sys_noise_est=Sys_noise_est_sum./Length_Names_train;



%% 
PSNR_sum=0;
SSIM_sum=0;
for num = 1 :Length_Names_test 
    image_name = File_test(num).name; 

    img_test=double(imread(strcat(test_path,image_name)));

    img=img_test(:,:,1:image_size_c)+sys_noise+5*randn(size( sys_noise ));

    img_rem_sysnoise=img-Sys_noise_est;

    [peaksnr, ssim] = MSIQA(img_test, img_rem_sysnoise);

    PSNR_sum=PSNR_sum+peaksnr;
    SSIM_sum=SSIM_sum+ssim;


end
PSNR_average=PSNR_sum/Length_Names_test;
SSIM_average=SSIM_sum/Length_Names_test;

