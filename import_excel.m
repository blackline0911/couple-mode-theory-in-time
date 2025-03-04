function varargout = import_excel(file_path,sheet, inspect_lambda_min, inspect_lambda_max)
%     inspect_lambda_min,  inspect_lambda_max should be nano meter    
    sheet_list = sheetnames(file_path);
    a = false;
    while (~a)
        for i = 1:numel(sheet_list)
            if (strcmp(sheet_list(i),sheet))
                a = true;
            end
        end
        if a~=true
            disp('the sheet name does not exist in the excel file');
            return;
        end
    end
    if (a)
        data = readmatrix(file_path,'Sheet',sheet);
        wavelength = data(:,1);
        L = length(wavelength);
        wl_min = min(wavelength);
        wl_max = max(wavelength);
        index_min = floor((inspect_lambda_min*1e-9-wl_min)/(wl_max-wl_min)*L);
        index_max = ceil((inspect_lambda_max*1e-9-wl_min)/(wl_max-wl_min)*L);
        for i=1:size(data,2)
            varargout{i} = data(index_min:index_max,i);
        end
    end
