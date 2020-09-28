fid1=fopen('1111.txt','r');
fid2=fopen('newtext.txt','w');
sign=1;
mm=1;

if fid1==-1
    error('The txt file is not available!');
else
    m=input('Input the total number of all the data:');
    n=input('Input the number of data which need to be modified:');

    id=unique(round(rand(1,n)*m));
    while length(id)<n
        id=unique([id,round(rand()*m)]);
    end
   
    frewind(fid1);
    while ~feof(fid1)    
        tline = fgets(fid1); % String of each line in the file.
        if length(tline)>6
            det=tline(1:6);
            if strcmp(det,'CQUAD4')
                tem=str2num(tline(10:15));
                if length(tem)==1 && tem>=100000 
                    if mm<=n && sign==id(mm)
                        fprintf('The modified line is: Line%d\n',id(mm));% Test for the line.
                        tline(18)='2';
                        mm=mm+1;
                    end            
                    sign=sign+1;
                end
                
            elseif length(tline)>45 && strcmp(tline(1:45),'$ Elements and Element Properties for region')
                break;
            end
            
        end
        fprintf(fid2,tline);        
    end
    
    
    
    while ~feof(fid1)
        tline = fgets(fid1);
        fprintf(fid2,tline); 
    end

    status1=fclose(fid1);
    status2=fclose(fid2);
    if status1==0 && status2==0
        disp('The file is closed.');
    else
        disp('The file is not closed');
    end

end

% x0=2;
% ite=10;
% 
% 
% f=@(x) (x*x+x-10);
% df=@(x) (2*x+1);
% 
% for i=1:ite
%     T=1;
%     miu=1;
%     if df(x0)==0
%         break;
%     end
%     x1=x0-miu*f(x0)/df(x0);
%     while abs(f(x1))>=abs(f(x0))
%         miu=miu/2;
%         x1=x0-miu*f(x0)/df(x0);
%         T=T+1;
%         if T>50
%             break;
%         end
%     end
%     x0=x1;
% end

