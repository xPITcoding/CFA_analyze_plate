function [o_mask,results]=regiongrowing(I,I_org,crit,cluster_size,threshold)

%
%I = 16bit greay value image
%crit = maximum allowed difference in grayvalues to the power of 2
%cluster_size = minimal amount of pixel to count as a cell cluster
%treshold = in orginal image

n_all = 0;
n_big_enough = 0;
n_round_enough = 0;

format compact;
mask=zeros(size(I));
dim=size(I);
o_mask=zeros(size(I));


I_org=im2double(I_org);
objid=0;
results=[];

h=waitbar(0,'analyzing well ...');

for rs=1:dim(1)
    waitbar(rs/dim(1),h);
    for cs=1:dim(2)        
        if I(rs,cs)==1 && mask(rs,cs)==0;
            objid=objid+1;
            point_count=1;
            pnt_vec=zeros(1,3);
            pnt_vec(1,:)=[0 rs cs];
            mask(rs,cs)=objid;
            val=I(rs,cs);
            A=[];
            n_all = n_all + 1;
            
            while point_count>0
                p=pnt_vec(1,2:3);
                A=[A;p];
                
                point_count=point_count-1;
                
                for r=p(1)-1:p(1)+1
                    for c=p(2)-1:p(2)+1
                        if (r>0 && c>0 && r<=dim(1) && c<=dim(2) && mask(r,c) == 0)
                            quaddiff=abs(val-I(r,c));
                            if (quaddiff<crit)
                                pnt_vec=[pnt_vec;[quaddiff r c]];
                                point_count=point_count+1;
                                
                                mask(r,c)=objid;
                            end
                        end
                    end
                end
                if (point_count>0)
                    pnt_vec=pnt_vec(2:point_count+1,:);
                end
            end
            
            %size(A,1)
            
            if size(A,1)<cluster_size
               for r=1:size(A,1)
                   mask(A(r,1),A(r,2))=0;
                   I(A(r,1),A(r,2))=0;
               end
               objid=objid-1;
            else
                n_big_enough = n_big_enough + 1;
                opnt_vec=[];
                for j=1:size(A,1)
                    v=I_org(A(j,1),A(j,2));
                    
                    if (v>threshold)
                        o_mask(A(j,1),A(j,2))=1;
                        opnt_vec=[opnt_vec;[v A(j,1) A(j,2)]];
                        
                    end
                end
                
                %so hier rundheitsdings was nichts macht
                [V,D]=eig(cov(opnt_vec(:,3),opnt_vec(:,2)));
                round = D(1,1)/D(2,2);
                [K,Area]=convhull(opnt_vec(:,3),opnt_vec(:,2));
                compactness = Area / length(opnt_vec);
                if ((abs(round-1)<0.75 | Area>5*cluster_size)  & compactness > 0.55) % 0.1 w�re weniger rund als 0.01 usw.
                    n_round_enough = n_round_enough + 1;
                    results(objid,1)=length(opnt_vec);
                    results(objid,2)=Area;
                    results(objid,3)=sum(opnt_vec(:,1))/length(opnt_vec);
                    results(objid,4)=sum(opnt_vec(:,1));
                    results(objid,5)=round; % habs jetzt auch noch als parameter angef�gt
                else 
                   % hier wird rausgeschmissen wenn es nicht rund ist
                   for r=1:size(A,1)
                       mask(A(r,1),A(r,2))=0;
                       o_mask(A(r,1),A(r,2))=0;
                       I(A(r,1),A(r,2))=0;
                   end
                   objid=objid-1;
                end
            end
            
            %size(results,1)
            
            % here we can now calculate statistics
            % N = number of pixel
            % A = area of the convex hull of the pixel 
            % B = average brightness of the pixel ... [K,B] = convhull(X,Y)
            % I = integrated brightness
                   
        end
    end
end

disp(sprintf('%d clusters found, %d larger than %d pixel, %d round shaped',n_all,n_big_enough,cluster_size,n_round_enough));

close(h);
            
           
                    
                    
                 
                            
                