clear all
clc

%% 
train_path = 'data\';
test_path = '';

File_train = dir(strcat(train_path,'*.tif')); 
Length_Names_train = length(File_train);
disp(Length_Names_train)


File_test = dir(strcat(test_path,'*.tif'));  
Length_Names_test = length(File_test);
disp(Length_Names_test)


%% 
image_size_h=*;
image_size_w=*;
image_size_c=*;

%% 
Sigma=65; %system noise level
sys_noise=Sigma*randn(image_size_h,image_size_w,image_size_c);
Sys_noise_est=zeros(size(sys_noise));
iter=zeros(1,image_size_c);
img_list=zeros(image_size_h,image_size_w,image_size_c,Length_Names_train);
for in = 1 :Length_Names_train 
    image_name = File_train(in).name;
    img_orig=double(imread(strcat(train_path,image_name)));
    img_list(:,:,:,in)=img_orig(:,:,1:image_size_c)+sys_noise+5*randn(size(sys_noise));
end
Opts=Opts_Set (Sigma);
for cc=1:image_size_c
    [Sys_noise_est(:,:,cc) , iter(cc)]     =     SC_GRSC_Sysnoise( squeeze(img_list(:,:,cc,:)),Opts);
end

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

