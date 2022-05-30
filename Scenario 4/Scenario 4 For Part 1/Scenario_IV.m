%% Zuhal Ceren Aydın  2444370
%% Burak Baki Bakır   2530103
clear all
clc

T = 20;
N = 240;



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
l_infs  = length(coords(1,240*0.05*0.5+1:240*0.05)) ;
l_is = length(coords(1,1:240*0.05*0.5));
infs = [coords(1,240*0.05*0.5+1:240*0.05) ; coords(2,240*0.05*0.5+1:240*0.05) ; zeros(1,l_infs)];
isos = [coords(1,1:240*0.05*0.5) ; coords(2,1:240*0.05*0.5) ; zeros(1,l_is)];
isos_i = isos(1:2,:);
not_infs = [coords(1,240*0.05+1:240);coords(2,240*0.05+1:240)];
iter_inf = zeros(1,length(infs(1,:)));
iter_n = zeros(1,length(not_infs(1,:)));
iter_isos = zeros(1,length(isos(1,:)));
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
vac_pi = [];
sec_vac = [];
n_v = [];
n_tv = [];
n_iv = [];
n_ndd = [];
n_dv = 0;
t_it = 1:150;
it_v = 20:150;
oy = true;
oh = true;
total_i = length(infs(1,:))+length(isos(1,:));
nn_v =0;
total_h = 0;
total_vi = 0;

for z = 1:150
    n_nd = 0;
    n_ni = 0;
    n_nh = 0;
    n_ive = 0;
    nn_v = 0;

    %%%%% vaccination process
     %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
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
        if oy == true
            vac_p = [vac_p(:,:) not_infs(:,1:floor(delta_3*ln_i))];
            vac_p(3:4,:) = [zeros(1,length(not_infs(2,1:floor(delta_3*ln_i)))) ; ones(1,length(not_infs(2,1:floor(delta_3*ln_i))))];
            oy = false;
        else
            vac_p = [vac_p(1,:) not_infs(1,1:floor(delta_3*ln_i)) ; vac_p(2,:) not_infs(2,1:floor(delta_3*ln_i)) ; vac_p(3,:) zeros(1,length(not_infs(2,1:floor(delta_3*ln_i)))) ; vac_p(4,:) ones(1,length(not_infs(2,1:floor(delta_3*ln_i))))];
        end
        nn_v = nn_v+ length(not_infs(:,1:floor(delta_3*ln_i)));
        n_v = [n_v  nn_v];
        
        not_infs(:,1:floor(delta_3*ln_i)) = [];

        t_vi = t_vi + nn_v;
        %%%%%%%%% walk vaccinated people
        k = 1;
        if length(vac_p(1,:)) >= 1 
            while k <= length(vac_p(1,:))
                coor = [vac_p(1,k)  vac_p(2,k)];
                new_c = randomWalk(coor);
                vac_p(1,k) = new_c(1);
                vac_p(2,k) = new_c(2);
                vac_p(3,k) = vac_p(3,k) +1;
                if vac_p(3,k) >= 3
                    if rand >= 0.2
                        sec_vac = [sec_vac(:,:) vac_p(1:2,k)];
                        vac_p(:,k) = [];
                        k = k-1;
                    else
                        vac_p(4,k) = 0;
                    end
                end
                k = k++1;
            end
        end
        n_tv = [n_tv  t_vi];
    end
    %%%%%%%%% walk nonisolated people
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

       %%%%%%%%% walk isolated people
    if length(isos(1,:)) >= 1 
        for k = 1:length(isos(1,:))
            coor = [isos_i(1,k)  isos_i(2,k)];
            new_c = isoWalk(coor);
            isos(1,k) = new_c(1);
            isos(2,k) = new_c(2);
            iter_isos(1,k) = iter_isos(1,k)+1;
        end
    end
    %%%%%%%% walk healed
   if isempty(heal) == false 
       if length(heal(1,:)) >= 1 
            for k = 1:length(heal(1,:))
                coor = [heal(1,k)  heal(2,k)];
                new_c = randomWalk(coor);
                heal(1,k) = new_c(1);
                heal(2,k) = new_c(2);
                heal(4,k) = heal(4,k);
                if heal(3,k) == 1 && heal(4,k) == 3
                    rnd = rand;
                    if rnd >= 0.2
                        sec_vac = [sec_vac(:,:) heal(1:2,k)];
                        heal(:,k) = [];
                    end
                end
            end
       end
   end

    %%%%%%%%% infection checker for vaccinated
     %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    if z >= 20
        k = 1;
        while k <= length(vac_p(1,:)) 
           tf_c = true;
           o = 1;
           ta_c = true;
           while o <= length(infs(1,:)) && ta_c == true
               if vac_p(1,k) == infs(1,o) && vac_p(2,k) == infs(2,o)
                   pr = rand;
                   pr2 = rand;
                   if vac_p(4,k) == 1
                       if pr <= 0.05 % get infected
                            n_ive = n_ive +1;
                            if pr2 >= 0.5
                                infs = [infs(1,:) vac_p(1,k) ; infs(2,:) vac_p(2,k) ; infs(3,:) 1];
                                iter_inf = [iter_inf 0];
                                n_ni =  n_ni +1;
                            else
                                isos = [isos(1,:) vac_p(1,k) ; isos(2,:) vac_p(2,k) ; isos(3,:) 1];
                                isos_i = [isos_i(1,:) vac_p(1,k) ; isos_i(2,:) vac_p(2,k)];
                                iter_isos = [iter_isos 0];
                                n_ni =  n_ni +1;
                            end                         
                           vac_p(:,k) = [];                   
                           k = k-1;
                           tf_c = false;
                       end
                   else
                       if pr <= 0.5 % get infected
                            n_ive = n_ive +1;
                            if pr2 >= 0.5
                                infs = [infs(1,:) vac_p(1,k) ; infs(2,:) vac_p(2,k) ; infs(3,:) 1];
                                iter_inf = [iter_inf 0];
                                n_ni =  n_ni +1;
                            else
                                isos = [isos(1,:) vac_p(1,k) ; isos(2,:) vac_p(2,k) ; isos(3,:) 1];
                                isos_i = [isos_i(1,:) vac_p(1,k) ; isos_i(2,:) vac_p(2,k)];
                                iter_isos = [iter_isos 0];
                                n_ni =  n_ni +1;
                            end                         
                           vac_p(:,k) = [];                   
                           k = k-1;
                           tf_c = false;
                       end
                   end
                   ta_c = false;
               end
               o = o+1;
           end

                %%%%%%%%% infection checker with isolated for vaccinated
                %%%%%%%%% people
           if tf_c == true
                o = 1;
                while o <= length(isos(1,:)) && ta_c == true 
                    if vac_p(1,k) == isos(1,o) && vac_p(2,k) == isos(2,o)
                        pr = rand;
                        pr2 = rand;
                        if vac_p(4,k) == 1
                            if pr <= 0.05 % get infected
                                n_ive = n_ive +1;
                                if pr2 >= 0.5
                                    infs = [infs(1,:) vac_p(1,k) ; infs(2,:) vac_p(2,k); infs(3,:) 1];
                                    iter_inf = [iter_inf 0];
                                    n_ni =  n_ni +1;
                                else
                                    isos = [isos(1,:) vac_p(1,k) ; isos(2,:) vac_p(2,k); isos(3,:) 1];
                                    isos_i = [isos_i(1,:) vac_p(1,k) ; isos_i(2,:) vac_p(2,k)];
                                    iter_isos = [iter_isos 0];
                                    n_ni =  n_ni +1;
                                end
                                vac_p(:,k) = [];
                                k = k-1;
                            end
                        else
                            if pr <= 0.5 % get infected
                                n_ive = n_ive +1;
                                if pr2 >= 0.5
                                    infs = [infs(1,:) vac_p(1,k) ; infs(2,:) vac_p(2,k); infs(3,:) 1];
                                    iter_inf = [iter_inf 0];
                                    n_ni =  n_ni +1;
                                else
                                    isos = [isos(1,:) vac_p(1,k) ; isos(2,:) vac_p(2,k); isos(3,:) 1];
                                    isos_i = [isos_i(1,:) vac_p(1,k) ; isos_i(2,:) vac_p(2,k)];
                                    iter_isos = [iter_isos 0];
                                    n_ni =  n_ni +1;
                                end
                                vac_p(:,k) = [];
                                k = k-1;
                            end
                        end
                        ta_c = false;
                    end
                    o = o+1;
                end
           end

           k = k+1;
        end
    end

    %%%%%%%%%%%%%%%% infection checker for nonvaccinated people
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    k=1;
    while k <= length(not_infs(1,:)) 
        tf_c = true;
        o = 1;
        ta_c = true;
        while o <= length(infs(1,:)) && ta_c == true
            if not_infs(1,k) == infs(1,o) && not_infs(2,k) == infs(2,o)
                pr = rand;
                pr2 = rand;
                if pr >= 0.5 % get infected
                    if pr2 >= 0.5
                        infs = [infs(1,:) not_infs(1,k) ; infs(2,:) not_infs(2,k) ; infs(3,:) 0];
                        iter_inf = [iter_inf 0];
                        n_ni =  n_ni +1;
                    else
                        isos = [isos(1,:) not_infs(1,k) ; isos(2,:) not_infs(2,k) ; isos(3,:) 0];
                        isos_i = [isos_i(1,:) not_infs(1,k) ; isos_i(2,:) not_infs(2,k)];
                        iter_isos = [iter_isos 0];
                        n_ni =  n_ni +1;
                    end
                    not_infs(:,k) = [];                   
                    k = k-1;
                    tf_c = false;
                end
                ta_c = false;
            end
            o = o+1;
        end
        %%%%%%%%% infection checker with isolated
        if tf_c == true
            o = 1;
            while o <= length(isos(1,:)) && ta_c == true
                if not_infs(1,k) == isos(1,o) && not_infs(2,k) == isos(2,o)
                    pr = rand;
                    pr2 = rand;
                    if pr >= 0.5 % get infected
                        if pr2 >= 0.5
                            infs = [infs(1,:) not_infs(1,k) ; infs(2,:) not_infs(2,k) ; infs(3,:) 0];
                            iter_inf = [iter_inf 0];
                            n_ni =  n_ni +1;
                        else
                            isos = [isos(1,:) not_infs(1,k) ; isos(2,:) not_infs(2,k) ; isos(3,:) 0];
                            isos_i = [isos_i(1,:) not_infs(1,k) ; isos_i(2,:) not_infs(2,k) ];
                            iter_isos = [iter_isos 0];
                            n_ni =  n_ni +1;
                        end
                        not_infs(:,k) = [];
                        k = k-1;
                    end
                    ta_c = false;
                end
                o = o+1;
            end
        end
        k = k+1;
    end

    %healing process
    j = 1;
    while j <= length(iter_inf) 
        if iter_inf(j)==30
            if oh == true
                        heal = zeros(4,1);
            end
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
                if infs(3,j) == 1  
                    heal = [heal(1,:) infs(1,j) ;heal(2,:) infs(2,j) ; heal(3,:) 1 ; heal(4,:) 0];
                else
                    heal = [heal(1,:) infs(1,j) ;heal(2,:) infs(2,j) ; heal(3,:) 0 ; heal(4,:) 0];
                end
                infs(:,j) = [];
                iter_inf(j) = []; 
                n_nh =  n_nh +1;
                j = j-1;
                if oh == true
                   heal(:,1) = [];
                   oh = false;
                end
            end
        end
        j = j+1;
    end



        %%%isolated healing process
    j = 1;
    while  j <= length(iter_isos) 
        if iter_isos(j)==30
            p_h = rand;
            if p_h > 0.95 % die
                if isos(3,j) == 1
                    n_dv = n_dv +1;
                end
                isos(:,j) = [];
                iter_isos(j) = [];
                isos_i(:,j) = [];
                n_d = n_d+1;
                j = j-1; 
                n_nd = n_nd +1;
            elseif p_h <= 0.95 %healed
                if isos(3,j) == 1
                    heal = [heal(1,:) isos(1,j) ;heal(2,:) isos(2,j) ; heal(3,:) 1 ; heal(4,:) 0];
                else
                    heal = [heal(1,:) isos(1,j) ;heal(2,:) isos(2,j) ; heal(3,:) 0 ; heal(4,:) 0];
                end
                isos(:,j) = [];
                iter_isos(j) = [];
                isos_i(:,j) = [];
                j = j-1;
                n_nh =  n_nh +1;
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
    

tiledlayout(5,2);


nexttile
plot(t_it,t_nni,"r",'LineWidth',2)
title('Number of newly infected people in each iteration')
xlim([0 175])
ylim([0 35])

nexttile
plot(t_it,t_ni,"r",'LineWidth',2)
title('Total number of infected people in the system in each iteration')
xlim([0 175])
ylim([0 240])


nexttile
plot(t_it,t_nnh,"g",'LineWidth',2)
title('Number of newly healed people in each iteration')
xlim([0 175])
ylim([0 30])

nexttile
plot(t_it,t_nh,"g",'LineWidth',2)
title('Total number of healed people in the system in each iteration')
xlim([0 175])
ylim([0 240])


nexttile
plot(t_it,t_nnd,'k','LineWidth',2)
title('Number of people that are died in that iteration')
xlim([0 175])
ylim([0 10])

nexttile
plot(t_it,t_nd,'k','LineWidth',2)
title('Total number of dead people in the system in each iteration')
xlim([0 175])
ylim([0 40])

nexttile
plot(it_v,n_v,'b','LineWidth',2)
title('Number of vaccinated people in each iteration,')
xlim([0 175])
ylim([0 120])

nexttile
plot(it_v,n_tv,'b','LineWidth',2)
title('Total number of vaccinated people in the system in each iteration')
xlim([0 175])
ylim([0 240])

nexttile
plot(t_it,n_iv,'r','LineWidth',2)
title('Number of infected people although they get vaccinated in each iteration,')
xlim([0 175])
ylim([0 20])

nexttile
plot(t_it,n_ndd,'k','LineWidth',2)
title('Number of dead people who got vaccinated')
ylabel('Number of dead people who got vaccinated')
xlabel('Number of iterations')
xlim([0 175])
ylim([0 8])
hold on




