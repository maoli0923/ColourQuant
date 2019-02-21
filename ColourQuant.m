%This is the code for an article "ColourQuant: a high-throughput technique to extract and quantify colour phenotypes from plant images" in a book chapter 

%3. Method
%3.1 Image Acquisition
%%%5.Colour correct images
%%%%%Balance the RGB color for each image so that white swatch R=255 G=255 B=255;

file=dir('*.jpg'); %if the format is png, then use '*.png' instead
for k=1:length(file)
    I=imread(file(k).name); %load the image
    %imtool(I) can be used to view the image and use mouse to point the white swatch region to know its position
    %Note that the (x,y) in the image refer to (w2,w1)
    w1=3440:3500;w2=554:742; 
    %This is the position of the white swatch
    %If this position differs image by image, extra algorithm to
    %detect it should to be added. For example, search the regions and find the region having the
    %colour closer to white. Or using imtool(I) mannually finds the region
    %image by image
    R=I(w1,w2,1);
    R=reshape(double(R),size(R,1)*size(R,2),1);
    G=I(w1,w2,2);
    G=reshape(double(G),size(G,1)*size(G,2),1);
    B=I(w1,w2,3);
    B=reshape(double(B),size(B,1)*size(B,2),1);
    I(:,:,1)=double(I(:,:,1))+255-floor(mean(R));
    I(:,:,2)=double(I(:,:,2))+255-floor(mean(G));
    I(:,:,3)=double(I(:,:,3))+255-floor(mean(B));
    imwrite(I,['C_' file(k).name(1:end-4) '.jpg']); %adding "C" in the filename to avoid replace the raw
end

%3.2. Object Segmentation
file=dir('C*.jpg'); %load the colour corrected images
for k=1:length(file)
    I=imread(file(k).name);
    %imtool(I) can be used to know the region of the color guide and other
    %tags
    II=imcrop(I,[1 270 3000 3000]);%crop out the color guide and other tags
    HSV=rgb2hsv(II);
    [a b]=find(HSV(:,:,2)>0.15); % the threshold value need to be customized for each experiment 
    BBW=zeros(size(II,1),size(II,2));
    for j=1:length(a)
        BBW(a(j),b(j))=1;
    end
    CC = bwconncomp(BBW);
    N=CC.PixelIdxList;
    Ind=[];
    for j=1:length(N)
        if length(N{j})>5000
            Ind=[Ind j];
        end
    end
    for j=1:length(Ind)
        III=zeros(size(BBW));
        III(N{Ind(j)})=1;
        [x y]=find(III==1);
        BW=imcrop(III,[min(y)-2,min(x)-2,max(y)-min(y)+4,max(x)-min(x)+4]);
        BW=imfill(BW,'holes'); %for each objest, this is the binary image
        h=imcrop(II,[min(y)-2,min(x)-2,max(y)-min(y)+4,max(x)-min(x)+4]);
        hR=h(:,:,1);
        hG=h(:,:,2);
        hB=h(:,:,3);
        w=find(BW==0); %background
        hR(w)=255;hG(w)=255;hB(w)=255;
        h=cat(3,hR,hG,hB);
        h=uint8(h);
        filename=sprintf(['Seg_' file(k).name(1:end-4),'_','%d','.mat'],j);
        save(filename, 'h')
    end
end

%3.3 Color analysis
%3.3.1 Mean and variance
file=dir('Seg_*.mat'); %load segmented image
Lab_Mean=zeros(length(file),3);
Lab_Var=zeros(length(file),3);
Lmin=[];Lmax=[];amin=[];amax=[];bmin=[];bmax=[]; %this is for section 3.3.2
for k=1:length(file)  
    load(file(k).name)
    lab=rgb2lab(h);%convert to lab color
    clear V
    V(:,1)=reshape(lab(:,:,1),size(h,1)*size(h,2),1);
    V(:,2)=reshape(lab(:,:,2),size(h,1)*size(h,2),1);
    V(:,3)=reshape(lab(:,:,3),size(h,1)*size(h,2),1);
    w=find(V(:,1)==100&V(:,2)==0&V(:,3)==0);%white background, it can also be computed using binary image
    V(w,:)=[];
    Lab_Mean(k,:)=mean(V);
    Lab_Var(k,:)=var(V);
    Lmin=[Lmin min(V(:,1))];Lmax=[Lmax max(V(:,1))];
    amin=[amin min(V(:,2))];amax=[amax max(V(:,2))];
    bmin=[bmin min(V(:,3))];bmax=[bmax max(V(:,3))];
end

%3.3.2 Gaussian density estimator
boundary=[min(Lmin) max(Lmax) min(amin) max(amax) min(bmin) max(bmax)];
xi=[floor(boundary(1))-10:4:ceil(boundary(2))+10];
yi=[floor(boundary(3))-10:4:ceil(boundary(4))+10];
zi=[floor(boundary(5))-10:4:ceil(boundary(6))+10];
bandwidth=[4 4 4]; %bandwith can be custormized, smaller will result higher resolution, but more computational time.
file=dir('Seg_*.mat'); %load segmented image
clear LabGD
for k=1:length(file)
    load(file(k).name)   
    lab=rgb2lab(h);%convert to lab color
    clear V
    V(:,1)=reshape(lab(:,:,1),size(h,1)*size(h,2),1);
    V(:,2)=reshape(lab(:,:,2),size(h,1)*size(h,2),1);
    V(:,3)=reshape(lab(:,:,3),size(h,1)*size(h,2),1);
    w=find(V(:,1)==100&V(:,2)==0&V(:,3)==0);%white background, it can also be computed using binary image
    V(w,:)=[];
    f=ksdensity3d(V,xi,yi,zi,bandwidth);
    LabGD(k,:)=reshape(f,1,length(xi)*length(yi)*length(zi));
end

%for different zones: border and center
file=dir('Seg_*.mat'); %load segmented image
clear BorderGD
clear CenterGD
for k=1:length(file)
    load(file(k).name)
    w=find(h(:,:,1)==255&h(:,:,2)==255&h(:,:,3)==255);
    BW=ones(size(h,1),size(h,2));
    BW(w)=0;
    [start_r,start_c] = find(BW,1,'first');
    V = bwtraceboundary(BW,[start_r start_c],'W',8,Inf,'counterclockwise');
    [x y]=find(BW==1);
    X=[x y];
    DD=pdist2(X,V);
    Dis=zeros(1,length(X));
    for jj=1:length(X)
        Dis(jj)=min(DD(jj,:));
    end
    VV=V;
    VV=VV-repmat(mean(VV),length(VV),1);
    Scale=norm(VV,'fro')/sqrt(length(VV));
    BorderInd=find(Dis<0.15*Scale); %This can be customized 
    CenterInd=find(Dis>max(Dis)-0.25*Scale); %This can be customized
    lab=rgb2lab(h);
    clear V
    clear U
    V(:,1)=reshape(lab(:,:,1),size(h,1)*size(h,2),1);
    V(:,2)=reshape(lab(:,:,2),size(h,1)*size(h,2),1);
    V(:,3)=reshape(lab(:,:,3),size(h,1)*size(h,2),1);
    w=find(BW==1);
    U=V(w(BorderInd),:);
    [f,xi,yi,zi]=ksdensity3d(U,xi,yi,zi,bandwidth);
    BorderGD(k,:)=reshape(f,1,length(xi)*length(yi)*length(zi));
    U=V(w(CenterInd),:);
    [f,xi,yi,zi]=ksdensity3d(U,xi,yi,zi,bandwidth);
    CenterGD(k,:)=reshape(f,1,length(xi)*length(yi)*length(zi));
end

%Distance based on different zones
dfull=dist(LabGD');
dborder=dist(BorderGD');
dcenter=dist(CenterGD');
D=sqrt(dfull.^2+dborder.^2+dcenter.^2);

%3.3.3 Circular deformation
file=dir('Seg_*.mat'); %load segmented image
lambda=0.00001; %a parameter for TPS
for k=1:length(file)
    load(file(k).name)
    landmarks=readPoints(h);  %mouse click on the tip and based and hit "Enter"
    w=find(h(:,:,1)==255&h(:,:,2)==255&h(:,:,3)==255);
    BW=ones(size(h,1),size(h,2));
    BW(w)=0;
    [start_r,start_c] = find(BW,1,'first');
    V = bwtraceboundary(BW,[start_r start_c],'W',8,Inf,'counterclockwise');
    theta=atan(-(landmarks(1,2)-mean(V(:,2)))/(landmarks(2,2)-mean(V(:,1))));
    Rot=[cos(theta) sin(theta);-sin(theta) cos(theta)]; %alignment
    [x y]=find(BW==1);
    score0=[x y];
    score0=score0-repmat(mean(V),length(score0),1);
    V=V-repmat(mean(V),length(V),1);
    V=V*Rot;
    score0=score0*Rot;
    score0=score0/norm(V,'fro')*sqrt(length(V));
    V=V/norm(V,'fro')*sqrt(length(V));
    [theta,r] = cart2pol(V(:,1),V(:,2));
    clear U
    [U(:,1) U(:,2)]=pol2cart(theta,1.5);
    ctrl0=V;
    ctrl1=U;
    [K] = tps_compute_dist_matrix(ctrl0);
    [Pm,Q1,Q2,R] = tps_set_matrices(ctrl0);
    [affine non_affine] = tps_compute_params(Q1,Q2,R,ctrl1,K,lambda);
    new_pt_temp=tps_defomed_pts(score0,affine,non_affine,ctrl0);
    VV=new_pt_temp(:,2:end);
    clear C
    C(:,1)=reshape(h(:,:,1),size(h,1)*size(h,2),1);
    C(:,2)=reshape(h(:,:,2),size(h,1)*size(h,2),1);
    C(:,3)=reshape(h(:,:,3),size(h,1)*size(h,2),1);
    w=find(C(:,1)==255&C(:,2)==255&C(:,3)==255);
    C(w,:)=[];
    scatter(VV(:,2),-VV(:,1),10,double(C)/255,'filled');axis equal off
    set(gca, 'Color', 'none');
    filename=sprintf(['Circular_' file(k).name(1:end-4) '.png']);
    export_fig(filename)% save to transparent .png
    close(gcf)
end

file=dir('Circular*.png');
for k=1:length(file)
    I=imread(file(k).name);
    I=imresize(I,[70 70]);
    lab=rgb2lab(I);
    clear V
    V(:,1)=reshape(lab(:,:,1),size(I,1)*size(I,2),1);
    V(:,2)=reshape(lab(:,:,2),size(I,1)*size(I,2),1);
    V(:,3)=reshape(lab(:,:,3),size(I,1)*size(I,2),1);
    X(k,:)=[V(:,1);V(:,2);V(:,3)]';
end

%perform PCA and get eigen leaves
[COEFF SCORE latent]=pca(X);
MeanX=mean(X);
for j=1:2
    PC1p=3*std(SCORE(:,j))*COEFF(:,j)'+MeanX;
    V=reshape(PC1p,size(X,2)/3,3);
    II(:,:,1)=reshape(V(:,1),size(I,1),size(I,2));
    II(:,:,2)=reshape(V(:,2),size(I,1),size(I,2));
    II(:,:,3)=reshape(V(:,3),size(I,1),size(I,2));
    RGB=lab2rgb(II);
    imtool(RGB)
    PC1p=-3*std(SCORE(:,j))*COEFF(:,j)'+MeanX;
    V=reshape(PC1p,size(X,2)/3,3);
    II(:,:,1)=reshape(V(:,1),size(I,1),size(I,2));
    II(:,:,2)=reshape(V(:,2),size(I,1),size(I,2));
    II(:,:,3)=reshape(V(:,3),size(I,1),size(I,2));
    RGB=lab2rgb(II);
    imtool(RGB)
end
    