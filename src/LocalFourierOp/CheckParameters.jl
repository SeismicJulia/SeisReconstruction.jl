
function check_parameters(in,Nw,Nkx,Nky,Noverlap_t,Noverlap_x, Noverlap_y;analysis_win=nothing)
    
    
    #I should use a loop to check the parameters at once and print the right infoarmation (associated to the right variable)
    #if something fails.
    #=
     if typeof(size(in)) != Tuple{Int64, Int64}
        message="The input must be at least 2D array!"
         error(message)
     end
    =# 
    
     if Nw == 0.0
         message="Nw must be greater than 1."
         error(message)
     else
        
        Nw=Int(Nw)
        if Nw < 1
            message="Nw must be greater than 1"
            error(message)
        end
     end
  
     if Noverlap_t > Nw
        message="Noverlap in time must be less than Nw."
        error(message)
        #Noverlap = Nw ÷ 2
     else
        Noverlap_t = Int(Noverlap_t)#toma el int
        if Noverlap_t >= Nw
            message=" Noverlap must be less than Nw. The condition Noverlap < Nw is not accomplished."
            error(message)
        end
     end

     hop_t= Nw- Noverlap_t;
    
    
      if Nkx == 0.0
         message="Nkx must be greater than 1."
         error(message)
     else
        
        Nkx=Int(Nkx)
        if Nkx < 1
            message="Nk must be greater than 1"
            error(message)
        end
     end
  
     if Noverlap_x > Nkx
        message="Noverlap in space must be less than Nk."
        error(message)
        #Noverlap = Nw ÷ 2
     else
        Noverlap_x = Int(Noverlap_x)#toma el int
        if Noverlap_x >= Nkx
            message=" Noverlap must be less than Nk. The condition Noverlap < Nk is not accomplished."
            error(message)
        end
     end

     hop_x= Nkx- Noverlap_x;
    
    
    #y
    
    if Nky == 0.0
         message="Nky must be greater than 1."
         error(message)
     else
        
        Nky=Int(Nky)
        if Nky < 1
            message="Nk must be greater than 1"
            error(message)
        end
     end
  
     if Noverlap_y > Nky
        message="Noverlap in space must be less than Nky."
        error(message)
        #Noverlap = Nw ÷ 2
     else
        Noverlap_y = Int(Noverlap_y)#toma el int
        if Noverlap_y >= Nky
            message=" Noverlap must be less than Nky. The condition Noverlap < Nk is not accomplished."
            error(message)
        end
     end

     hop_y= Nky- Noverlap_y;
    
    
    
    return Int(hop_t),Int(hop_x), Int(hop_y)

end