%[hutest,hltest] = envelope(Y);
%x_sim = linspace(0,1-(1/100000),100000);
%Y_sim = 2.187e4*sin(2*pi*225*x_sim - pi/3.1);%-1.411e4;
%Y_sim_neg = 2.187e4*sin(2*pi*225*x_sim + pi/1.45);%-1.411e4;

peak = [];
peak_norm = [];
T_temp = [];
temp = [];
peak_pos = [];
peak_neg = [];
bvmean_neg = 0;
bvmean_pos = 0;


  %for g = 1:size(dV_sig_orig,2) 
      
  %      h = 1;
   
  temp = mean(dV_sig_orig,2);
  BV_u_mean = mean(BV_u);
  BV_l_mean = mean(BV_l);
  
        for i = 2:size(dV_sig_orig,1)-1
        
            if temp(i) > 0 && temp(i) > temp(i-1) && temp(i) > temp(i+1)
            
                T_temp(h) = T(i);    
            
                peak_norm(h) = temp(i)/BV_u_mean;
            
                h = h + 1;
                
            end

        if temp(i) < 0 && temp(i) < temp(i-1) && temp(i) < temp(i+1)
            
            T_temp(h) = T(i);    
            
            peak_norm(h) = temp(i)/BV_l_mean;
            
            h = h + 1;
            
        end
        
    end

  
%  for i = 1:length(T_temp)
      
%      if peak_norm(i,g) > 0 && T_temp(i) > -20
          
%          temp(i) = 
%          bvmean_pos = 
          
          
  
  %{
  for m = 1:length(peak_norm,2)
  for o = 1:length(peak_norm,1)-1
        
        if T_temp(o+1,m) - T_temp(o,m) < 3
            
            T_temp(o,m) = mean([T_temp(o,m) T_temp(o+1,m)]);
            peak_norm(o,m) = mean([peak_norm(o,m) peak_norm(o+1,m)]);
            
        end
        
  end
  end
  %}
