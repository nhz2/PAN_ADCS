function init=sensor_prop(init)
%adds sensor parameters to init

init.IMU0_DCM_BODY= eye(3);
%DCM to translate from body basis to IMU0 basis
init.IMU0_MAG_FS= 4.0;
%The fullscale of the magnetometer is +- this (Gauss)
init.IMU0_MAG_GN= 6842.0;
%what is read when 1 gauss is applied
init.IMU0_MAG_ODR= 80.0;
%Output data rate in Hz also used for the temp sensor
init.IMU0_ZGAUSS= zeros([3 1]);
%gauss applied to get 0 output
init.IMU0_RMS= transpose([3.2e-3,3.2e-3,4.1e-3]);
%RMS noise in gauss
