f = fopen('u.data');
A=fscanf(f,'%d%d%d%d');
A=vec2mat(A,4);
fclose(f);
train=zeros(max(A(:,1)),max(A(:,2)));
test=zeros(max(A(:,1)),max(A(:,2)));
user_count=zeros(943);
ucnt=0;
cnt=0;
testcnt=0;


for x=A'                        %this loop is for creating 200 users train, test sets- fast- takes ~30 minutes
    if user_count(x(1))==0
        if ucnt>=200
            continue
        end
        ucnt=ucnt+1;
    end
    cnt=cnt+1;
    if user_count(x(1))<50 || testcnt>=5600
        train(x(1),x(2))=x(3);
        user_count(x(1))=user_count(x(1))+1;
    else
        testcnt=testcnt+1;
        test(x(1),x(2))=x(3);
    end
end

% for x=A'                            %this loop is for creating 943 users train, test sets- takes tiem- ~1 hour
%     cnt=cnt+1;
%     if user_count(x(1))<50 || testcnt>=20000
%         train(x(1),x(2))=x(3);
%         user_count(x(1))=user_count(x(1))+1;
%     else
%         testcnt=testcnt+1;
%         test(x(1),x(2))=x(3);
%     end
% end


testcnt
cnt