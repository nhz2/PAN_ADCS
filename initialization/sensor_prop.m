function init=sensor_prop(init)
%adds sensor parameters to init

init.IMU0_DCM_BODY= eye(3);
%DCM to translate from body basis to IMU0 basis
init.IMU0_POS_BODY= zeros([3 1]);
%position of IMU0 in m in body coordinates
init.IMU0_MAG_FS= 4.0;
%The fullscale of the magnetometer is +- this (Gauss)
init.IMU0_MAG_GN= 6842.0;
%what is read when 1 gauss is applied
init.IMU0_MAG_ODR= 10.0;
%Mag output data rate in Hz also used for the temp sensor
init.IMU0_MAG_ZGAUSS= zeros([3 1]);
%gauss applied to get 0 output
init.IMU0_MAG_RMS= transpose([3.2e-3,3.2e-3,4.1e-3]);
%RMS noise in gauss
init.IMU0_G_ODR= 100.0;
%Gyro output data rate in Hz
init.IMU0_G_TYOFF= zeros([3 1]);%transpose([0.17,0.17,0.17]);
%gyro offset, what the gyro reads in rad/s when not rotating
init.IMU0_G_GN= 13096.0;
%what is read when 1 rad/s
init.IMU0_G_FS= 2.182;
%The fullscale of the gyro is +- this (rad/s)
init.IMU0_G_RMS= transpose([1.2e-3,1.2e-3,1.2e-3]);
%RMS noise in rad/s