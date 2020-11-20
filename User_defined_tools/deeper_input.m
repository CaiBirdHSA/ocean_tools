function deeper_input(input_data)

temp=ncread(input_data,'temp');
salt=ncread(input_data,'salt');
[im,jm,kb,nt]=size(temp);

mask3d=1-isnan(temp(:,:,:,1));  %land is 0, ocean is 1;
flag = 0;

for i=1:im
    for j=1:jm
        for k=1:kb-1
            if(mask3d(i,j,k)==1 && mask3d(i,j,k+1)==0 )
                kstart= k+1;
                kend  = min(k+6,kb);
                temp(i,j,kstart:kend,:) = repmat(temp(i,j,k,:),1,1,kend-kstart+1,1);
                salt(i,j,kstart:kend,:) = repmat(salt(i,j,k,:),1,1,kend-kstart+1,1);
                flag=flag+1;
            end
        end
    end
end
                    
if flag>0
    disp(['update the input_data: ', input_data]);
    ncwrite(input_data,'temp',temp);
    ncwrite(input_data,'salt',salt);
end

end