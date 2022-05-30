
%% Zuhal Ceren Aydın  2444370
%% Burak Baki Bakır   2530103
clear all
clc

T = 20;
N = 240;
avt_nni = zeros(200,150);
avt_ni = zeros(200,150);
avt_nnh = zeros(200,150);
avt_nh = zeros(200,150);
avt_nnd = zeros(200,150);
avt_nd = zeros(200,150);
n_ddd = zeros(200,150);
avt_ndd= zeros(200,150);



for uq = 1:200
    
    
    coords = zeros(2,240);
    infs = [];
    grid = zeros(T,T);
    
    for i= 1:240
        pl = true;
        while pl == true
            [a ,b] = coord();
            na = false;
            for k = 1:length(coords)
                if a == coords(1,k) && b == coords(2,k)
                    na= true;
                end
            end
            if na == false    
                coords(1,i) =  a;
                coords(2,i) =  b;
                pl=false;
            end
        end
    end
    l_infs= 240*0.05;
    infs = [coords(1,1:l_infs) ; coords(2,1:l_infs) ; zeros(1,l_infs)];
    not_infs = [coords(1,240*0.05+1:240);coords(2,240*0.05+1:240)];
    iter_inf = zeros(1,length(infs(1,:)));
    iter_n = zeros(1,length(not_infs(1,:)));
    heal = [];
    n_d = 0;
    t_vi = 0;
    t_nni = [];
    t_ni = [];
    t_nnh = [];
    t_nh = [];
    t_nnd = [];
    t_nd = [];
    vac_p = [];
    n_v = [];
    n_tv = [];
    n_iv = [];
    n_ndd = [];
    n_dv = 0;
    
    t_it = 1:150;
    it_v = 20:150;
    total_i = length(infs(1,:));
    total_h = 0;
    total_vi = 0;
    
    for z = 1:150
        n_nd = 0;
        n_ni = 0;
        n_nh = 0;
        n_ive = 0;
        nn_v= 0;
    
        %%%%% vaccination process
        if z >= 20
            delta_3 = 1/(2*(z-19));
            ln_i = length(not_infs(1,:));
            ln_h = 0;
            if isempty(heal) == false 
                for urt = 1:length(heal(1,:))
                    if heal(3,urt) == 0
                        ln_h = ln_h+1;
                    end
                end
                ln_h = floor(ln_h*delta_3);
                nn_v = nn_v +ln_h;
                urt = 1;
                while urt<=ln_h
                    if heal(3,urt) == 0
                        heal(3,urt) = 1;
                        urt = urt+1;
                    else
                        ln_h = ln_h+1;
                        urt = urt+1;
                    end 
                end
            end
            nn_v = nn_v+ length(not_infs(:,1:floor(delta_3*ln_i)));
            n_v = [n_v  nn_v];
           
            vac_p = [vac_p(:,:) not_infs(:,1:floor(delta_3*ln_i))];
            not_infs(:,1:floor(delta_3*ln_i)) = [];
            t_vi = t_vi + nn_v;
            %%%%%%%%% walk vaccinated people
            if length(vac_p(1,:)) >= 1 
                for k = 1:length(vac_p(1,:))
                    coor = [vac_p(1,k)  vac_p(2,k)];
                    new_c = randomWalk(coor);
                    vac_p(1,k) = new_c(1);
                    vac_p(2,k) = new_c(2);
                end
            end
            
    
            n_tv = [n_tv  t_vi];
        end
        %%%%%%%%% walk normal  people
        if length(infs(1,:)) >= 1 
            for k = 1:length(infs(1,:))
                coor = [infs(1,k)  infs(2,k)];
                new_c = randomWalk(coor);
                infs(1,k) = new_c(1);
                infs(2,k) = new_c(2);
                iter_inf(1,k) = iter_inf(1,k)+1;
            end
        end
        if length(not_infs(1,:)) >= 1 
            for k = 1:length(not_infs(1,:))
                 coor = [not_infs(1,k)  not_infs(2,k)];
                 new_c = randomWalk(coor);
                 not_infs(1,k) = new_c(1);
                 not_infs(2,k) = new_c(2);
                 iter_n(1,k) = iter_n(1,k) +1;
            end
        end
    
        %%%%%%%%% infection checker for vaccinated
        if z >= 20
            k = 1;
            while k <= length(vac_p(1,:)) 
               o = 1;
               ta_c = true;
               while o <= length(infs(1,:)) && ta_c == true
                   if vac_p(1,k) == infs(1,o) && vac_p(2,k) == infs(2,o)
                       pr = rand;
                       if pr <= 0.05 % get infected
                           infs = [infs(1,:) vac_p(1,k) ; infs(2,:) vac_p(2,k) ; infs(3,:) 1];
                           iter_inf = [iter_inf 0];
                           n_ni =  n_ni +1;
                           vac_p(:,k) = [];   
                           n_ive = n_ive +1;
                           k = k-1;
                       end
                       ta_c = false;
                   end
                   o = o+1;
               end
               k = k+1;
            end
            
        end
        %%%%%%%%% infection checker for nonvaccinated
        k=1;
        while k <= length(not_infs(1,:)) 
            o = 1;
            ta_c = true;
            while o <= length(infs(1,:)) && ta_c == true
                if not_infs(1,k) == infs(1,o) && not_infs(2,k) == infs(2,o)
                    pr = rand;
                    if pr >= 0.5 % get infected
                        infs = [infs(1,:) not_infs(1,k) ; infs(2,:) not_infs(2,k) ; infs(3,:) 0];
                        iter_inf = [iter_inf 0];
                        iter_n(:, k) = [];
                        n_ni =  n_ni +1;
                        not_infs(:,k) = [];                   
                        k = k-1;
                    end
                    ta_c = false;
                end
                o = o+1;
            end
            k = k+1;
            
        end
        %healing process
        j = 1;
        while j <= length(iter_inf) 
            if iter_inf(j)==30
                p_h = rand;
                if p_h > 0.95 % die
                    if infs(3,j) == 1
                        n_dv = n_dv +1;
                    end
                    infs(:,j) = [];
                    iter_inf(j) = [];
                    n_d = n_d+1;
                    j = j-1; 
                    n_nd = n_nd +1;
                elseif p_h <= 0.95 % heal
                    heal = [heal(:,:) infs(:,j) ];
                    infs(:,j) = [];
                    iter_inf(j) = []; 
                    n_nh =  n_nh +1;
                    j = j-1;
                end
            end
            j = j+1;
        end
        %%%datas needed
        total_i = total_i +n_ni;
        t_ni = [t_ni total_i];
        total_h = total_h +n_nh;
        t_nh = [t_nh total_h];
        total_vi = n_ive+total_vi;
        n_iv = [n_iv total_vi];
        t_nd = [t_nd n_d];
        t_nni = [t_nni n_ni];
        t_nnh = [t_nnh n_nh];
        t_nnd = [t_nnd n_nd];
        n_ndd = [n_ndd n_dv];
    
        
        end
    avt_ndd(uq,:) = n_ndd;
    avt_nni(uq,:) = t_nni;
    avt_ni(uq,:) = t_ni;
    avt_nnh(uq,:) = t_nnh;
    avt_nh(uq,:) = t_nh;
    avt_nnd(uq,:) = t_nnd;
    avt_nd(uq,:) = t_nd;
    n_ativ(uq,:) = n_iv;
    n_atv(uq,:) = n_v;
    n_attv(uq,:) = n_tv;
    
        
end
avtt_nni = zeros(1,150);
avtt_ni = zeros(1,150);
avtt_nnh = zeros(1,150);
avtt_nh = zeros(1,150);
avtt_nnd = zeros(1,150);
avtt_nd = zeros(1,150);
n_aativ = zeros(1,150);
n_aatv = zeros(1,131);
n_aattv = zeros(1,131);
n_nnndd = zeros(1,150);
for zt = 1:150
    avtt_nni(1,zt) = mean(avt_nni(:,zt));
    avtt_ni(1,zt) = mean(avt_ni(:,zt));
    avtt_nnh(1,zt) = mean(avt_nnh(:,zt));
    avtt_nh(1,zt) = mean(avt_nh(:,zt));
    avtt_nnd(1,zt) = mean(avt_nnd(:,zt));
    avtt_nd(1,zt) = mean(avt_nd(:,zt));
    n_aativ(1,zt) = mean(n_ativ(:,zt));
    n_nnndd(1,zt) = mean(avt_ndd(:,zt));
end
for ztt = 1:131
    n_aatv(1,ztt) = mean(n_atv(:,ztt));
    n_aattv(1,ztt) = mean(n_attv(:,ztt));
end



figure(2)
plot(t_it,avtt_nni,"m",'LineWidth',2)
title('Number of newly infected people in each iteration')
ylabel('Number of newly infected people')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 35])
hold on

figure(3)
plot(t_it,avtt_ni,"m",'LineWidth',2)
title('Total number of infected people in the system in each iteration')
ylabel('Total number of infected people in the system')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 240])
hold on


figure(4)
plot(t_it,avtt_nnh,"m",'LineWidth',2)
title('Number of newly healed people in each iteration')
ylabel('Number of newly healed people')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 30])
hold on

figure(5)
plot(t_it,avtt_nh,"m",'LineWidth',2)
title('Total number of healed people in the system in each iteration')
ylabel('Total number of healed people in the system')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 240])

hold on


figure(6)
plot(t_it,avtt_nnd,'m','LineWidth',2)
title('Number of people that are died in that iteration')
ylabel('Number of people that are died')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 10])
hold on

figure(7)
plot(t_it,avtt_nd,'m','LineWidth',2)
title('Total number of dead people in the system in each iteration')
ylabel('Total number of dead people in the system')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 40])
hold on

figure(8)
plot(it_v,n_aatv,"m",'LineWidth',2)
title('Number of vaccinated people in each iteration,')
ylabel('Number of vaccinated people')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 120])
hold on

figure(9)
plot(it_v,n_aattv,'m','LineWidth',2)
title('Total number of vaccinated people in the system in each iteration')
ylabel('Total number of vaccinated people in the system')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 240])
hold on

figure(10)
plot(t_it,n_aativ,'m','LineWidth',2)
title('Number of infected people although they get vaccinated in each iteration,')
ylabel('Number of infected people although they get vaccinated')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 20])
hold on

figure(11)
plot(t_it,n_nnndd,'m','LineWidth',2)
title('Number of dead people who got vaccinated')
ylabel('Number of dead people who got vaccinated')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 8])
hold on
