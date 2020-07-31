
curFig             = scf(100001);
clf(curFig,"reset");
//demo_viewCode("membrane.sce");

drawlater();

xselect(); //raise the graphic window
ntt=2;
nt=2409;
i=2;

// set a new colormap
//-------------------
cmap= curFig.color_map; //preserve old setting
curFig.color_map = jetcolormap(64);

outfile=sprintf('out/output%d.out',i);

u(:,:,ntt)=fscanfMat(outfile, "%lg");

u(1,1,2)=-1;
u(1,2,2)=1;
plot3d1(x,y,u(:,:,2),35,45,' ',ebox=[0 0.8 0 0.1 -2.0 2.0]);
u(1,1,2)=0;
u(1,2,2)=0;
s=gce(); //the handle on the surface
s.color_flag=1 ; //assign facet color according to Z value
s.color_mode=-1;
title("evolution of a 3d surface","fontsize",3)



drawnow();



for i=2000:nt
  outfile=sprintf('out/output%d.out',i);
  u(:,:,ntt)=fscanfMat(outfile, "%lg");
  
  
  realtime(i); //wait till date 0.1*i seconds
  //s.data.z = (sin((I(i)/10)*x)'*cos((I(i)/10)*y))';
  s.data.z = u(:,:,ntt);
  //s.data_bounds=[0,0,-2;1,1,2];
end

