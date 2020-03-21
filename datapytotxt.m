load('pythondata.mat')
lgn = 0.90 * reshape(gl',1,69120*150);
save('gl.txt','lgn','-ASCII');
