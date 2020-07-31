%
% Scilab ( http:%www.scilab.org/ ) - This file is part of Scilab
% Copyright (C) INRIA
%
% This file is distributed under the same license as the Scilab package.
%

% =============================================================================
% Demonstrate animation based on the evolution of a 3D surface
% ============================================================================




% set a new colormap
%-------------------
%cmap= curFig.color_map; %preserve old setting
%curFig.color_map = jetcolormap(64);

%The initial surface definition 
%----------------------
x=linspace(-pi,pi,50);
y=x;
Z=15*sin(x)'*cos(y);
myones=ones(50,50);
[mmx,mmy]=meshgrid(x,y);
%Creates and set graphical entities which represent the surface
%--------------------------------------------------------------
plot3d1(x,y,Z,35,45,' ');
s=gce(); %the handle on the surface
s.color_flag=1 ; %assign facet color according to Z value
title("evolution of a 3d surface","fontsize",3)

h=surf(Z);                              %Create initial surface return a handle for this graphic
axis([0 50 0 50 -10 10]);               %define axis and keep fixed 
hold on;                                %hold the graphics



I=4000:-0.1:1;
%realtimeinit(0.1);;%set time step (0.1 seconds)  and date reference


%drawnow();

for i=1:max(size(I))
  %realtime(i); %wait till date 0.1*i seconds
  %s.data.z = (sin((I(i)/10)*x)'*cos((I(i)/10)*y))';
  %s.data.z = sin(6*pi*i/max(size(I)))*(sin((2)*x)'*sin((4)*y))';
  %s.data.z = sin((100*9*pi*i*myones/max(size(I))-2*mmx+4*mmy));%+sin((100*8*pi*i*myones/max(size(I))-4*mmx+2*mmy))+sin((100*20*pi*i*myones/max(size(I))-2*mmx+20*mmy));
  Z = 5.0*sin((50*4.5*pi*i*myones/max(size(I)))).*(sin((mmx/2)).*sin((mmy/2)))+4.0*sin((100*4.5*pi*i*myones/max(size(I)))).*(sin((mmx)).*sin((mmy)))+2.0*sin((100*9*pi*i*myones/max(size(I)))).*(sin((-2*mmx)).*sin((+4*mmy)))+4.0*sin((100*45*pi*i*myones/max(size(I)))).*(sin((3*mmx)).*sin((-3*mmy)))+0.5*sin((100*80*pi*i*myones/max(size(I)))).*(sin((6*mmx)).*sin((-6*mmy)))+0.25*sin((100*160*pi*i*myones/max(size(I)))).*(sin((12*mmx)).*sin((-12*mmy)))+0.2*sin((100*180*pi*i*myones/max(size(I)))).*(sin((24*mmx)).*sin((-24*mmy)))+3*sin((100*180*pi*i*myones/max(size(I)))).*(sin((6*mmx)).*sin((-5*mmy)))+3*sin((100*360*pi*i*myones/max(size(I)))).*(sin((-12*mmx)).*sin((10*mmy)));;


    
    pause(0.15);
    set(h,'ZData',Z)                    %reset the z data for the grahics window with grahics handle h
    drawnow

end
