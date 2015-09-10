function predictions = movie_prediction(train, test)

test_users=any(test')';
L=size(test_users);
L=L(1);
M=size(train);
M=M(2);
resultmat=zeros(size(test));

sth=0.35;

avg_r=zeros(L);
cnt=1;
for x=train'              %computing average user ratings
    cnt2=0;
    sum=0;
    for y=x'
        if y==0
            continue;
        end
        sum=sum+y;
        cnt2=cnt2+1;
    end
    if cnt2~=0
        avg_r(cnt)=sum/cnt2;
    end
    cnt=cnt+1;
end


S=zeros(L,L);
Si=zeros(M,M);
'computing user similarity'
for i=1:L
    for j=1:i
        sum1=0;
        sum2=0;
        sum3=0;
        for k=1:M
            if train(i,k)==0 || train(j,k)==0
                continue;
            end
            sum1=sum1+(train(i,k)-avg_r(i))*(train(j,k)-avg_r(j));
            sum2=sum2+(train(i,k)-avg_r(i))*(train(i,k)-avg_r(i));
            sum3=sum3+(train(j,k)-avg_r(j))*(train(j,k)-avg_r(j));
        end
        if sum2==0 || sum3==0
            S(i,j)=0;
            S(j,i)=0;
        else
            S(i,j)=sum1/(sqrt(sum2)*sqrt(sum3));
            S(j,i)=Si(i,j); 
        end
    end
end
'computing item similarity'
for i=1:M
    for j=1:i
        sum1=0;
        sum2=0;
        sum3=0;
        for k=1:L
            if train(k,i)==0 || train(k,j)==0
                continue;
            end
            sum1=sum1+(train(k,i)-avg_r(k))*(train(k,j)-avg_r(k));
            sum2=sum2+(train(k,i)-avg_r(k))*(train(k,i)-avg_r(k));
            sum3=sum3+(train(k,j)-avg_r(k))*(train(k,j)-avg_r(k));
        end
        if sum2==0 || sum3==0
            Si(i,j)=0;
            Si(j,i)=0;
        else
            Si(i,j)=sum1/(sqrt(sum2)*sqrt(sum3));
            Si(j,i)=Si(i,j); 
        end
    end
end

A=zeros(M,M);
for i=1:M
    [sorted,I]=sort(Si(i,:),'descend');
    for j=1:5
        if sorted(1,j)>sth
            A(i,I(j))=1;
        end
    end
    
end    

true_cnt=0;
total_cnt=0;
ucnt=0;

mae=0;
rmse=0;

for i=1:L                            %iterating over active users
    if test_users(i)~=1
        continue
    end
    ucnt=ucnt+1;

    G=A;
    for j=1:M
        if test(i,j)==0
            G(j,:)=0;
            G(:,j)=0;
        end
    end
    
    y=zeros(M);
    for j=1:M                  %computing values for local evidence nodes y_i
        if test(i,j)==0
            continue;
        end
        y(j)=avg_r(i);
        sum1=0;
        sum2=0;
        
        for k=1:L
            if train(k,j)==0
                continue;
            end
            sum1=sum1+S(i,k)*(train(k,j)-avg_r(k));
            sum2=sum2+abs(S(i,k));
        end
        if sum2~=0 
            y(j)=y(j)+sum1/sum2;
        end
    end
    
    
    phi=zeros(M,5);
    for j=1:M                    %computing values for initial evidence phi_i
        if test(i,j)==0
            continue
        end
        if y(j)<1
            phi(j,1)=1;
        elseif y(j)>5
            phi(j,5)=1;
        elseif isinteger(y(j))
            phi(j,y(j))=1;
        else
            x=(y(j)-ceil(y(j)))/(floor(y(j))-ceil(y(j)));
            phi(j,floor(y(j)))=x;
            phi(j,ceil(y(j)))=1-x;
        end    
    end
    
    psi=zeros(M,M,25);
    for j=1:M                 %computing values for compatibility function or clique of size 2
        if test(i,j)==0
            continue
        end
        for k=1:M
            if test(i,k)==0
                continue
            end
            sigma=((1-Si(j,k))/(1-sth) + 1)/sqrt(2);
            sum=0;
            for ij=1:5
                for ik=1:5
                    psi(j,k,5*(ij-1)+ik)=exp(-((ij-ik)^2)/(sigma^2));
                    sum=sum+psi(j,k,5*(ij-1)+ik);
                end
            end
            psi(j,k,:)=psi(j,k,:)/sum;
        end
    end
    
    m=ones(M,M,5)*(1/5);
    cnt=0;
    p=zeros(M,5);
    result=zeros(M);
    
    while true   %belief propagation loop
        temp=ones(M,M,5)*(1/5);
        for a=1:M
            if test(i,a)==0
                continue
            end
            neighbours=find(G(a,:));
            for b=1:M
                if test(i,b)==0
                    continue
                end
                sum=0;
                for rb=1:5
                    sum2=0;
                    for ra=1:5
                        prod=1;
                        for n=neighbours
                            if n==b
                                continue
                            end
                            prod=prod*m(n,a,ra);
                        end
                        prod=prod*psi(a,b,5*(ra-1)+rb)*phi(a,ra);
                        sum2=sum2+prod;
                    end
                    temp(a,b,rb)=sum2;
                    sum=sum+sum2;
                end
                
                temp(a,b,:)=temp(a,b,:)/sum;
                if any(isnan(temp(a,b,:)))
                    i
                    a
                    b
                end
            end
        end
        
        if all((m-temp)<0.001)
            break;
        end
        if cnt>100
            break;
        end
        m=temp;
        cnt=cnt+1;
    end
    ucnt
    for j=1:M           %computing marginals and evaluation matrices
        if test(i,j)==0
            continue
        end
        neighbours=find(G(j,:));
        sum=0;
        for rj=1:5
            prod=phi(j,rj);
            for n=neighbours
                prod=prod*m(n,j,rj);
            end
            p(j,rj)=prod;
            sum=sum+prod;
        end
        p(j,:)=p(j,:)/sum;
        
        
        for rj=1:5
            result(j)=result(j)+rj*p(j,rj);
        end
        
        if isnan(result(j))
            continue;
        end
        
        resultmat(i,j)=result(j);
        
        if test(i,j)==round(result(j))
            true_cnt=true_cnt+1;
        end
        mae=mae+abs(test(i,j)-result(j));
        rmse=rmse+(test(i,j)-result(j))*(test(i,j)-result(j));
        total_cnt=total_cnt+1;
    end
    
end

mae
rmse

total_cnt

mae=mae/total_cnt
rmse=sqrt(rmse/total_cnt)


predictions=true_cnt/total_cnt

assignin('base','result',resultmat);
