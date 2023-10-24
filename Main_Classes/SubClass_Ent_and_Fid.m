classdef SubClass_Ent_and_Fid
%--------------------------------------------------------------------------
%Created by Eva Takou
%Last modified: Oct 24, 2023
%--------------------------------------------------------------------------

%This class has the properties: one tangles, Kraus, Fidelity
%entangling power of unitary (one-tangles, or M-way entanglement).
   
  properties
      
      Uval     {mustBeNumeric}  
      Ek            %Kraus operators (could be empty)
      Fid           %Fidelity (as calculated from the total system using the Kraus operators)
      maxBoundPi    %Max bound of one-tangles for CR evolutions
      one_tangles   %One-tangles for a general U (or CR-type)
      Max_Bound     %Max bound of one-tangles for general U.
      epM_uni       %M-way entangling power (of target subspace)
      epM_nonuni    %M-way entangling power (of target subspace)
      epM_dephased  %M-way entangling power (of target subspace)
      
  end
  
  
  properties (Constant, Hidden)
       
        Q=1/sqrt(2)*[1  0  0  1i ; 
             0 1i  1  0  ;
             0 1i -1  0  ;
             1  0  0  -1i];
  
  end
    
  methods (Static)  
         
      
      function Ek=Get_Kraus(U,Nnuc,Target_Spins)
          
          %Target_Spins: Subsystem indices for nuclear spins to keep.
          
          if any(Target_Spins>Nnuc)
              
              error('The nuclear spin indices exceed the # of total nuclei.')
              
          end
          
          K    = length(Target_Spins); %# of subsystems in target nuclear subspace
          Kall = K+1;                  %# of susbystems in target subspace (including electron)
          
          %------- Parameters for the environment -------------------------
          
          Nenv      = Nnuc-K;                   %nuclei in environment
          dim_E     = 2^Nenv;                   %dimension of environment
          id_E      = eye(dim_E);               %Identity
          ket0_E    = zeros(dim_E,1);          
          ket0_E(1) = 1;                        %initial state of environment
          e0_E      = kron(eye(2^Kall),ket0_E); %Promote to total space
          
          %--- Bring all target spins to first positions ------------------
          
          Target_Spins = sort(Target_Spins)+1; %Add 1 because 1st subsystem is electron
         
          if any(Target_Spins>Kall) %Need to permute

              for jj=2:length(Target_Spins)+1

                  if Target_Spins(jj-1)~=jj

                      SWAP = ArbDimSWAP(jj,Target_Spins(jj-1),Nnuc+1); 
                      U    = SWAP*U*SWAP'; %Permute to jj-th position

                  end

              end              
              
          end
          
          %"Trace-out" all the final subsystems
          
          Ek = zeros(2^(Kall),2^(Kall),dim_E);
          
          for k=1:dim_E
              
              ek = kron(eye(2^(Kall)),id_E(k,:)); %For systems we do not trace out, we need to put Id.
              Ek(:,:,k)=ek*U*e0_E;
              
          end
          
          
          %Check completeness of Kraus:
          compl=0;
          
          for k1=1:dim_E
              
              
              compl = compl + Ek(:,:,k1)'*Ek(:,:,k1);
                  
              
          end
          
          imag_compl = imag(compl);
          
          if all(imag_compl<1e-9)
              
              compl=real(compl);
              
          end
          
          if norm(compl-eye(2^(Kall)))>1e-9
              
              error('Completeness is not satisfied')
              
          end
          
          
      end
      
      function Infid=Get_Infid_From_Kraus(U0,Ek,Target_Spins)
          
          Kall = length(Target_Spins)+1;
          
          if log2(length(U0))~=Kall
              
             error('Dimension of target gate does not match the # of target susbystems.') 
             
          end
          
          expr1=0;
          expr2=0;
          
          for k=1:size(Ek,3)
             
              Mk=U0'*Ek(:,:,k);
              
              expr1=expr1+Mk'*Mk;
              expr2=expr2+abs(trace(Mk))^2;
              
              
          end
          
          m=2^Kall;
          
          Infid = 1 - 1/(m*(m+1))*(trace(expr1)+expr2);
          
          
          
          
          
      end
        
    function [Out]=Gate_Infid(Nnuc,Ntarget,phi0,phi1,n0,n1)
    %Input: Nnuc: total # of nuclear spins
    %       Ntarget: # of nuclear spins of target subspace
    %       phi0,phi1: Rot angles of unwanted nuclei
    %       n0,n1: Rot axes of unwanted nuclei
    %Output: Gate infidelity.

    %Analytical evaluation of gate infidelity which expresses deviation
    %of target gate (of target subsoace) from ideal due to the presence 
    %of unwanted nuclei.

     dTarget = 2^(Ntarget+1);   
     Nunw    = Nnuc-Ntarget;

     if Nunw~=length(phi0)

         error('Inconsistent # of unwanted spins and rot angles argument.')

     end

    %NumberOfKet ranges from 0 till 2^{L-K}-1 where L#of all nuclei
    %                                               K#of target ones

    %loop through all these combinations, read the m-th bit
    %and if the m-th position is 0, we are supposed to get
    %the factor cos(phi_0/2)-1i*nz0*sin(phi_0/2) for when the electron
    %is in 0 (and similarly for 1).
    %if the m-th position is 1 then we get the other coefficient

    expr=0;

    for NumberOfKet = 0:2^(Nunw)-1 

        CurrentKet=dec2bin(NumberOfKet,Nunw)-'0';

        %To get the logical indices do this:

        Nuc_in0 = CurrentKet==0;
        Nuc_in1 = CurrentKet==1;

        ang0_el0    = phi0(Nuc_in0); ang0_el1    = phi1(Nuc_in0);
        ang1_el0    = phi0(Nuc_in1); ang1_el1    = phi1(Nuc_in1);

        %for nuclei in |0>
        f0_ThisKraus_el0 = prod(cos(ang0_el0/2)-1i*(n0(Nuc_in0,3).').*sin(ang0_el0/2));
        f0_ThisKraus_el1 = prod(cos(ang0_el1/2)-1i*(n1(Nuc_in0,3).').*sin(ang0_el1/2));

        %for nuclei in |1>
        f1_ThisKraus_el0 = prod((-1i)*sin(ang1_el0/2).*(n0(Nuc_in1,1)+1i*n0(Nuc_in1,2)).');
        f1_ThisKraus_el1 = prod((-1i)*sin(ang1_el1/2).*(n1(Nuc_in1,1)+1i*n1(Nuc_in1,2)).');


        if  isempty(ang0_el0) && ~isempty(ang1_el0) %all nuclei in |1>

            f0_ThisKraus_el0=1;
            f0_ThisKraus_el1=1;

        elseif ~isempty(ang0_el0) && isempty(ang1_el0) %all nuclei in |0>

            f1_ThisKraus_el0=1;
            f1_ThisKraus_el1=1;

        end

         expr = expr + abs(f0_ThisKraus_el0*f1_ThisKraus_el0+...
                           f0_ThisKraus_el1*f1_ThisKraus_el1)^2;

    end


    F          = 1/(dTarget+1)*(1+2^(Ntarget-1)*expr); 
    Out.Infid  = 1-F;   

    end
        
        
  end
  
  methods  %M-way entanglement metrics

        function obj=NonUniMwayEp_CR(obj,Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1)
        
            %Input: Nnuc: Total # of nuclei
            %       Ntarget: # of target nuclei
            %       phi0,phi1: Rot angles of unwanted spins
            %       n0,n1: Rot axes of unwanted spins
            %       th0,th1: Rot angles of target spins
            %       u0,u1: Rot axes of target spins

            %To calculate the non-unitary M way entangling power of CR type
            %evolution operator. The non-unitarity arises due to partial trace
            %of unwanted nuclei. The metric includes the impact of unwanted
            %correlations on target M-way entangling power.
        
            if Ntarget~=length(th0)
           
                error('The length of rot angle array passed for target nuclei is not equal to # of target nuclei.') 
                
            end
            
            g1=@(indx) (cos(th0(indx)/2)*cos(th1(indx)/2)+...
                     dot(u0(indx,:),u1(indx,:))*sin(th0(indx)/2)*sin(th1(indx)/2))^2;
         
            d         = 2;
            prefactor = 1;        
         
            for ii=1:length(th0)
            
                prefactor = prefactor* (1-g1(ii));
            
            end
            
            prefactor = 1/2*prefactor *(d/(d+1))^(Ntarget+1);
            Nunw      = Nnuc-Ntarget;
        
            %NumberOfKet ranges from 0 till 2^{L-K}-1 where L#of all nuclei
            %                                               K#of target ones
        
            %Loop through all these combinations, read the m-th bit
            %and if the m-th position is 0, we are supposed to get
            %the factor cos(phi_0/2)-1i*nz0*sin(phi_0/2) for when the electron
            %is in 0 (and similarly for 1).
            %if the m-th position is 1 then we get the other coefficient
        
            %Do this for any combinations of 2 Kraus operators

            %coefficient of Kraus when the electron is in 0/1 and nuclei are in 0
            F0_S0 = @(PHIj,Spinsj)  prod(cos(PHIj/2)-1i*(n0(Spinsj,3).').*sin(PHIj/2));
            F1_S0 = @(PHIj,Spinsj)  prod(cos(PHIj/2)-1i*(n1(Spinsj,3).').*sin(PHIj/2));
        
            %coefficient of Kraus when the electron is in 0/1 and nuclei are in 1             
            F0_S1 = @(PHIj,Spinsj)  prod((-1i)*sin(PHIj/2).*(n0(Spinsj,1)+1i*n0(Spinsj,2)).');
            F1_S1 = @(PHIj,Spinsj)  prod((-1i)*sin(PHIj/2).*(n1(Spinsj,1)+1i*n1(Spinsj,2)).');
        
            expr  = 0;
            %extra = 0;
            
            parfor NumberOfKet = 0:2^(Nunw)-1 

                %disp(['Iter=',num2str(NumberOfKet),' out of:',num2str(2^(Nunw)-1)])
                
                CurrentKet=dec2bin(NumberOfKet,Nunw)-'0';
                Nuc_in0 = CurrentKet==0;  Nuc_in1 = CurrentKet==1;

                %for nuclei in |0>
                f0_ThisKraus_el0 = F0_S0(phi0(Nuc_in0),Nuc_in0);
                f0_ThisKraus_el1 = F1_S0(phi1(Nuc_in0),Nuc_in0);

                %for nuclei in |1>
                f1_ThisKraus_el0 = F0_S1(phi0(Nuc_in1),Nuc_in1);
                f1_ThisKraus_el1 = F1_S1(phi1(Nuc_in1),Nuc_in1);


                if  all(~Nuc_in0) && all(Nuc_in1) %all nuclei in |1>

                    f0_ThisKraus_el0=1;
                    f0_ThisKraus_el1=1;

                elseif all(Nuc_in0) && all(~Nuc_in1) %all nuclei in |0>

                    f1_ThisKraus_el0=1;
                    f1_ThisKraus_el1=1;

                end

                F0_ThisKraus = f0_ThisKraus_el0*f1_ThisKraus_el0;
                F1_ThisKraus = f0_ThisKraus_el1*f1_ThisKraus_el1;

                for NumberOfOtherKet = 0:2^(Nunw)-1 

                    OtherKet=dec2bin(NumberOfOtherKet,Nunw)-'0';

                    OtherNuc_in0 = OtherKet==0;  OtherNuc_in1 = OtherKet==1;

                    %for nuclei in |0>
                    f0_OtherKraus_el0 = F0_S0(phi0(OtherNuc_in0),OtherNuc_in0);
                    f0_OtherKraus_el1 = F1_S0(phi1(OtherNuc_in0),OtherNuc_in0);

                    %for nuclei in |1>
                    f1_OtherKraus_el0 = F0_S1(phi0(OtherNuc_in1),OtherNuc_in1);
                    f1_OtherKraus_el1 = F1_S1(phi1(OtherNuc_in1),OtherNuc_in1);

                    if  all(~OtherNuc_in0) && all(OtherNuc_in1) %all nuclei in |1>

                        f0_OtherKraus_el0=1;
                        f0_OtherKraus_el1=1;

                    elseif all(OtherNuc_in0) && all(~OtherNuc_in1) %all nuclei in |0>

                        f1_OtherKraus_el0=1;
                        f1_OtherKraus_el1=1;

                    end

                    F0_OtherKraus = f0_OtherKraus_el0*f1_OtherKraus_el0;
                    F1_OtherKraus = f0_OtherKraus_el1*f1_OtherKraus_el1;


                    %expr = expr + abs(F0_ThisKraus*F1_OtherKraus)^2 ...
                    %             + abs(F1_ThisKraus*F0_OtherKraus)^2;


                    expr = expr+real(conj(F0_ThisKraus)*conj(F1_OtherKraus)*...
                                     F1_ThisKraus*F0_OtherKraus);
                    
                    %extra = extra+2*real(F0_ThisKraus'*F1_OtherKraus'*F1_ThisKraus*F0_OtherKraus);


                end
            end
            
            expr = expr+1;
            expr = expr*prefactor;
            
        
            %if (-1)^(Ntarget+1)==1 %even target subspace.. extra factor of 1/2
                
            %    expr = expr+extra;
            %    expr=expr/2;

            %end
            
            %expr = prefactor*expr;
            obj.epM_nonuni=expr;   
       
        end

        function obj=NonUniMwayEp_CR_Faster(obj,Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1)
            %Input: Nnuc: Total # of nuclei
            %       Ntarget: # of target nuclei
            %       phi0,phi1: Rot angles of unwanted spins
            %       n0,n1: Rot axes of unwanted spins
            %       th0,th1: Rot angles of target spins
            %       u0,u1: Rot axes of target spins

            %To calculate the non-unitary M way entangling power of CR type
            %evolution operator. The non-unitarity arises due to partial trace
            %of unwanted nuclei. The metric includes the impact of unwanted
            %correlations on target M-way entangling power.
        
            if Ntarget~=length(th0)
           
                error('The length of rot angle array passed for target nuclei is not equal to # of target nuclei.') 
                
            end
            
            g1=@(indx) (cos(th0(indx)/2)*cos(th1(indx)/2)+...
                     dot(u0(indx,:),u1(indx,:))*sin(th0(indx)/2)*sin(th1(indx)/2))^2;
         
            d         = 2;
            prefactor = 1;        
         
            for ii=1:length(th0)
            
                prefactor = prefactor* (1-g1(ii));
            
            end
            
            prefactor = 1/2*prefactor *(d/(d+1))^(Ntarget+1);
            Nunw      = Nnuc-Ntarget;
        
            %NumberOfKet ranges from 0 till 2^{L-K}-1 where L#of all nuclei
            %                                               K#of target ones
        
            %Loop through all these combinations, read the m-th bit
            %and if the m-th position is 0, we are supposed to get
            %the factor cos(phi_0/2)-1i*nz0*sin(phi_0/2) for when the electron
            %is in 0 (and similarly for 1).
            %if the m-th position is 1 then we get the other coefficient
        
            %Do this for any combinations of 2 Kraus operators

            %coefficient of Kraus when the electron is in 0/1 and nuclei are in 0
            F0_S0 = @(PHIj,Spinsj)  prod(cos(PHIj/2)-1i*(n0(Spinsj,3).').*sin(PHIj/2));
            F1_S0 = @(PHIj,Spinsj)  prod(cos(PHIj/2)-1i*(n1(Spinsj,3).').*sin(PHIj/2));
        
            %coefficient of Kraus when the electron is in 0/1 and nuclei are in 1             
            F0_S1 = @(PHIj,Spinsj)  prod((-1i)*sin(PHIj/2).*(n0(Spinsj,1)+1i*n0(Spinsj,2)).');
            F1_S1 = @(PHIj,Spinsj)  prod((-1i)*sin(PHIj/2).*(n1(Spinsj,1)+1i*n1(Spinsj,2)).');
        
            expr  = 0;
            
            for NumberOfKet = 0:2^(Nunw)-1 

                disp(['Iter=',num2str(NumberOfKet),' out of:',num2str(2^(Nunw)-1)])
                
                CurrentKet=dec2bin(NumberOfKet,Nunw)-'0';
                Nuc_in0 = CurrentKet==0;  Nuc_in1 = CurrentKet==1;

                %for nuclei in |0>
                f0_ThisKraus_el0 = F0_S0(phi0(Nuc_in0),Nuc_in0);
                f0_ThisKraus_el1 = F1_S0(phi1(Nuc_in0),Nuc_in0);

                %for nuclei in |1>
                f1_ThisKraus_el0 = F0_S1(phi0(Nuc_in1),Nuc_in1);
                f1_ThisKraus_el1 = F1_S1(phi1(Nuc_in1),Nuc_in1);


                if  all(~Nuc_in0) && all(Nuc_in1) %all nuclei in |1>

                    f0_ThisKraus_el0=1;
                    f0_ThisKraus_el1=1;

                elseif all(Nuc_in0) && all(~Nuc_in1) %all nuclei in |0>

                    f1_ThisKraus_el0=1;
                    f1_ThisKraus_el1=1;

                end

                F0_ThisKraus = f0_ThisKraus_el0*f1_ThisKraus_el0;
                F1_ThisKraus = f0_ThisKraus_el1*f1_ThisKraus_el1;

                for NumberOfOtherKet = NumberOfKet+1:2^(Nunw)-1 

                    OtherKet=dec2bin(NumberOfOtherKet,Nunw)-'0';

                    OtherNuc_in0 = OtherKet==0;  OtherNuc_in1 = OtherKet==1;

                    %for nuclei in |0>
                    f0_OtherKraus_el0 = F0_S0(phi0(OtherNuc_in0),OtherNuc_in0);
                    f0_OtherKraus_el1 = F1_S0(phi1(OtherNuc_in0),OtherNuc_in0);

                    %for nuclei in |1>
                    f1_OtherKraus_el0 = F0_S1(phi0(OtherNuc_in1),OtherNuc_in1);
                    f1_OtherKraus_el1 = F1_S1(phi1(OtherNuc_in1),OtherNuc_in1);

                    if  all(~OtherNuc_in0) && all(OtherNuc_in1) %all nuclei in |1>

                        f0_OtherKraus_el0=1;
                        f0_OtherKraus_el1=1;

                    elseif all(OtherNuc_in0) && all(~OtherNuc_in1) %all nuclei in |0>

                        f1_OtherKraus_el0=1;
                        f1_OtherKraus_el1=1;

                    end

                    F0_OtherKraus = f0_OtherKraus_el0*f1_OtherKraus_el0;
                    F1_OtherKraus = f0_OtherKraus_el1*f1_OtherKraus_el1;


                    %expr = expr + abs(F0_ThisKraus*F1_OtherKraus)^2 ...
                    %             + abs(F1_ThisKraus*F0_OtherKraus)^2;


                    expr = expr+2*real(conj(F0_ThisKraus)*conj(F1_OtherKraus)*...
                                            F1_ThisKraus*F0_OtherKraus);
                    
                    

                end
                
                expr = expr + abs(F0_ThisKraus)^2*abs(F1_ThisKraus)^2;
                
            end
            
            expr = expr+1;
            expr = expr*prefactor;
        
            obj.epM_nonuni=expr;   
       
        end

        function obj=NonUniMwayEp_CR_Analytical(obj,Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1)
            %Input: Nnuc: Total # of nuclei
            %       Ntarget: # of target nuclei
            %       phi0,phi1: Rot angles of unwanted spins
            %       n0,n1: Rot axes of unwanted spins
            %       th0,th1: Rot angles of target spins
            %       u0,u1: Rot axes of target spins

            %To calculate the non-unitary M way entangling power of CR type
            %evolution operator. The non-unitarity arises due to partial trace
            %of unwanted nuclei. The metric includes the impact of unwanted
            %correlations on target M-way entangling power.
        
            if Ntarget~=length(th0)
           
                error('The length of rot angle array passed for target nuclei is not equal to # of target nuclei.') 
                
            end
            
            g1T=@(indx) (cos(th0(indx)/2)*cos(th1(indx)/2)+...
                     dot(u0(indx,:),u1(indx,:))*sin(th0(indx)/2)*sin(th1(indx)/2))^2;
         
            g1U=@(indx) (cos(phi0(indx)/2)*cos(phi1(indx)/2)+...
                     dot(n0(indx,:),n1(indx,:))*sin(phi0(indx)/2)*sin(phi1(indx)/2))^2;
                 
                 
            Nunw      = Nnuc-Ntarget;
        
            exprU=1;

            for jj=1:Nunw
                
                
                exprU=exprU*(g1U(jj) + (...
                      n0(jj,3)*sin(phi0(jj)/2)*cos(phi1(jj)/2)...
                     -n1(jj,3)*sin(phi1(jj)/2)*cos(phi0(jj)/2)...
                     -(n1(jj,1)*n0(jj,2)-n0(jj,1)*n1(jj,2))*sin(phi0(jj)/2)*sin(phi1(jj)/2)...
                             ...
                                 )^2 ...
                             );
                
            end
            
            exprU=(exprU+1)/2;
            
            %Calculate the target subspace entangling power:
            exprT=1;
            
            for jj=1:Ntarget
                
                exprT = exprT*(1-g1T(jj));
                
            end
            
            d         = 2;
            prefactor = (d/(d+1))^(Ntarget+1);
           
            obj.epM_nonuni=prefactor*exprT*exprU;   
       
        end
        
        function obj=NonUniMwayEp_CR_Aproximate(obj,Nnuc,Ntarget,phi0,phi1,n0,n1,th0,th1,u0,u1)
            %Input: Nnuc: Total # of nuclei
            %       Ntarget: # of target nuclei
            %       phi0,phi1: Rot angles of unwanted spins
            %       n0,n1: Rot axes of unwanted spins
            %       th0,th1: Rot angles of target spins
            %       u0,u1: Rot axes of target spins

            %To calculate the non-unitary M way entangling power of CR type
            %evolution operator. The non-unitarity arises due to partial trace
            %of unwanted nuclei. The metric includes the impact of unwanted
            %correlations on target M-way entangling power.
        
            if Ntarget~=length(th0)
           
                error('The length of rot angle array passed for target nuclei is not equal to # of target nuclei.') 
                
            end
            
            g1=@(indx) (cos(th0(indx)/2)*cos(th1(indx)/2)+...
                     dot(u0(indx,:),u1(indx,:))*sin(th0(indx)/2)*sin(th1(indx)/2))^2;
         
            d         = 2;
            prefactor = 1;        
         
            for ii=1:length(th0)
            
                prefactor = prefactor* (1-g1(ii));
            
            end
            
            g1_unwanted=@(indx) (cos(phi0(indx)/2)*cos(phi1(indx)/2)+...
                     dot(n0(indx,:),n1(indx,:))*sin(phi0(indx)/2)*sin(phi1(indx)/2))^2;
           
            expr=1;     
                 
            for jj=1:length(phi0)
                
                expr=expr*g1_unwanted(jj);
                
            end
                 
            prefactor = 1/2*prefactor *(d/(d+1))^(Ntarget+1);
            
            expr = prefactor*(expr+1);
            obj.epM_nonuni = expr;   
       
        end
        
        function obj=UniMwayEp_CR(obj,Ntarget,th0,th1,u0,u1)
        %Input:  Ntarget: # of target nuclei
        %        th0,th1: Rot angles of target nuclei
        %        u0,u1: Rot axes of target nuclei (Ntarget x 3 arrays)
        %Output: Unitary M-way entangling power of U, which is CR-type.
        
        d    = 2;
        M    = (Ntarget+1);
        pref = (d/(d+1))^M;
        
        G1_Makhlin=@(th0,th1,u0,u1) (cos(th0/2)*cos(th1/2)+dot(u0,u1)*sin(th0/2)*sin(th1/2))^2;
        
        expr = zeros(1,Ntarget);
        
        for jj=1:Ntarget
           
            expr(jj) = 1-G1_Makhlin(th0(jj),th1(jj),u0(jj,:),u1(jj,:));
            
        end
        
        obj.epM_uni = pref*prod(expr);
        
        
        end
        
        function obj=MWayEp(obj)
        
            %Output: M-way entangling power of arbitrary U that conists of even
            %# of qubits.
            
            U=obj.Uval;
        
            if isempty(U)
               error('Provide value of U on the object.') 
            end
        
            M = log2(size(U,1));
            d = 2;
            
            if (-1)^M~=1
                
                error('Can only evaluate for M even for non-CR evolutions.')
            end
            
            Om0p = 1;
            Pm   = 1;
            
            Pplus  = cell(1,M);
            Pminus = cell(1,M);
            
            for jj=1:M
               
                SWAPij     = ArbDimSWAP(jj,jj+M,2*M);
                Pplus{jj}  = 1/2*(eye(2^(2*M))+SWAPij);
                Pminus{jj} = 1/2*(eye(2^(2*M))-SWAPij);
                Om0p       = Om0p*Pplus{jj};
                Pm         = Pm*Pminus{jj};
                
            end
            
            Om0p        = Om0p/(d+1)^M;
            epM         = 2^M*trace(kron(U,U)*Om0p*kron(U',U')*Pm);
            obj.epM_uni = epM;
        end
        
        function obj=dephased_epM_CR(obj,n0,n1,phi0,phi1,lambda0,lambda1,theta,t,N)
        
            %Input:  n0,n1: Rot axes of target nuclei (passed as a list)
            %        phi0,phi1: Rot angles of target nuclei (passed as 1x
            %        Ntarget array)
            %        lambda0,lambda1: Dephasing probabilities
            %        theta: Dephasing angle (in units of inverse time)
            %        t: Time(s) of CPMG (or XY2) sequences
            %        N: Number of iterations (of CPMG sequences)
            %Output: Mway entangling power of CR-type U, consisting of control
            %qubit + target spins. The control expreriences dephasing, described by
            %the parameters lambda_j and theta (parameters of Kraus operators).
        
            M    = length(phi0)+1;
            pref = (2/3)^M;
            l0   = lambda0;
            l1   = lambda1;

            G1 = zeros(1,length(phi0));

            for jj=1:length(phi0)

                G1(jj) = (cos(phi0(jj)/2)*cos(phi1(jj)/2) + dot(n0{jj},n1{jj})* sin(phi0(jj)/2)*sin(phi1(jj)/2))^2;

            end
            
            G1prod = prod(1-G1);

        
            if length(N)==1 %We are not composing sequences of different unit times and iterations

                R = (l0^2+l1^2+2*l0*l1*cos(theta*t))^(2*N)  * (l0^2+l1^2 + 2*l0*l1*cos(2*theta*t))^N;
                epM = (1+R)/2*pref*G1prod;      

            elseif length(N)>1 %We are composing sequences of different unit times and iterations
                        
               R = 1;

               for jj=1:length(N)
                   
                  R = R *  ( l0^2 + l1^2 + 2 * l0 * l1 * cos(theta*t(jj)) )^(2*N(jj))  ...
                              *...
                           ( l0^2 + l1^2 + 2 * l0 * l1 * cos(2*theta*t(jj)) )^(N(jj));

               end

                epM = (1+R)/2*pref*G1prod;      

            end
            
            obj.epM_dephased = epM;
            
        end
        
        function obj=dephased_epM_brute_force(obj,wL,A,B,s0,s1,lambda0,lambda1,theta,t,N)
            %Function to calculate the M-tangling power 
            %with additional electronic dephasing.
            %Can only be applied to a CPMG sequence.
            
            l0 = lambda0;          
            l1 = lambda1;
            d  = 2;
            M    = length(A)+1;     %All M qubits
            Nnuc = M-1;             %# of nuclear spins
            Id   = eye(2^Nnuc);     %Identity matrix
            X    = [0 1 ; 1 0];
            XI   = kron(X,Id);
            Upi  = expm(-1i*pi/2*XI);
            
            %Electronic dephasing operators:
            K0=@(t)   sqrt(l0)*kron(diag([exp(1i*theta*t),exp(-1i*theta*t)]),Id);
            K1=@(t)   sqrt(l1)*kron(diag([exp(-1i*theta*t),exp(1i*theta*t)]),Id);
            
            temp           = SuperClass_Sequences(wL,A,B,s0,s1,Nnuc,1,1);
            [Hdiag,T,Tinv] = temp.Get_Ham;
            H0             = T*Hdiag*Tinv;
            
            %We now take all possible multiplications of the Kraus
            %operators. %2^{3*N} total # of Kraus operators.
            
            Kraus    = cell(1,2^(3*sum(N)));
            bitCombs = dec2bin(0:2^(3*sum(N))-1) - '0';
            
            for jj=1:size(bitCombs,1)
                
                TEMP = 1;
                cnt  = 1;
                T    = t(cnt);
                
                for kk=1:sum(N)

                    if kk==sum(N(1:cnt))+1 && length(t)>1
                        
                        cnt = cnt+1;    
                        T   = t(cnt);
                        
                    end
                    
                    m1   = bitCombs(jj,1+3*(kk-1));
                    m2   = bitCombs(jj,2+3*(kk-1));
                    m3   = bitCombs(jj,3+3*(kk-1));
                    K    = this_Kraus(m1,m2,m3,T,K0,K1,Upi,H0);
                    
                    TEMP = K*TEMP;
                    
                end
                
                Kraus{jj}=TEMP;
                
            end
            
            function K=this_Kraus(m1,m2,m3,t,K0,K1,Upi,H0)
                
                Uf = @(t,K) K*expm(-1i*H0*t);
                
                if m1==0 && m2==0 && m3==0
                    
                    K=Uf(t/4,K0(t/4))*Upi*Uf(t/2,K0(t/2))*Upi*Uf(t/4,K0(t/4));
                    
                elseif m1==0 && m2==0 && m3==1
                    
                    K=Uf(t/4,K0(t/4))*Upi*Uf(t/2,K0(t/2))*Upi*Uf(t/4,K1(t/4));

                elseif m1==0 && m2==1 && m3==0
                    
                    K=Uf(t/4,K0(t/4))*Upi*Uf(t/2,K1(t/2))*Upi*Uf(t/4,K0(t/4));
                    
                elseif m1==0 && m2==1 && m3==1

                    K=Uf(t/4,K0(t/4))*Upi*Uf(t/2,K1(t/2))*Upi*Uf(t/4,K1(t/4));
                    
                elseif m1==1 && m2==0 && m3==0
                    
                    K=Uf(t/4,K1(t/4))*Upi*Uf(t/2,K0(t/2))*Upi*Uf(t/4,K0(t/4));
                    
                elseif m1==1 && m2==0 && m3==1
                    
                    K=Uf(t/4,K1(t/4))*Upi*Uf(t/2,K0(t/2))*Upi*Uf(t/4,K1(t/4));

                elseif m1==1 && m2==1 && m3==0
                    
                    K=Uf(t/4,K1(t/4))*Upi*Uf(t/2,K1(t/2))*Upi*Uf(t/4,K0(t/4));
                    
                elseif m1==1 && m2==1 && m3==1
                    
                    K=Uf(t/4,K1(t/4))*Upi*Uf(t/2,K1(t/2))*Upi*Uf(t/4,K1(t/4));
                    
                end
                
                
            end
            
            Omp0=1;
            Pp  = cell(1,M);
            Pm  = cell(1,M);
            
            for jj=1:M
               
                Pp{jj} = 1/2*(eye(2^(2*M))+ArbDimSWAP(jj,jj+M,2*M));
                Pm{jj} = 1/2*(eye(2^(2*M))-ArbDimSWAP(jj,jj+M,2*M));
                Omp0   = Pp{jj}*Omp0;
                
            end
            
            
            if (-1)^M==1
                
                P=1;
                
                for jj=1:M
                   
                    P = P*Pm{jj};
                    
                end
                
            else
                
                P=Pp{1};
                
                for jj=2:M
                   
                    P = P*Pm{jj};
                    
                end                
                
            end
            
            Omp0 = Omp0/(d+1)^M;
            epM  = 0;

            for j1=1:length(Kraus)

                for j2=1:length(Kraus)

                    epM=epM+trace(kron(Kraus{j1},Kraus{j2})  * Omp0 * ...
                                  kron(Kraus{j1},Kraus{j2})' * P);

                end
            end            
            
            obj.epM_dephased = real(d^M* epM);
            
        end
        
  end
  
  methods %One-tangles
      
      function obj=one_tangles_CR(obj,Nnuc,phi0,phi1,n0,n1)
      %Calculate one-tangles for CR-type evolution operator.
      %Input: Nnuc: # of nuclear spins
      %       phi0,phi1: Rot angles
      %       n0,n1: Rot axes
      
          D  = Nnuc+1;
          G1 = @(phi0,phi1,n0n1) (cos(phi0/2)*cos(phi1/2)+n0n1*sin(phi0/2)*sin(phi1/2))^2;
          Fs = zeros(1,Nnuc);
          
          for ii=1:Nnuc
          
              Fs(ii)=G1(phi0(ii),phi1(ii),dot(n0(jj,:),n1(jj,1)));
        
          end
          
          epNuc                 = 2/9*(1-Fs);   
          ep_electron           = 1/3 - 1/3^(D)*prod(1+2*Fs);
          obj.one_tangles       = {epNuc,ep_electron};
          
          epmax_Nuc      = 2/9;
          epmax_Electron = 1/3-1/(3^D);
        
          obj.maxBoundPi = {epmax_Nuc,epmax_Electron,[]};
          
      end
      
      function obj=one_tangles_General_U_Zanardi(obj)
      %To calculate one-tangles of arbitrary U. We follow a method based on
      %Zanardi.
      
          U     = obj.Uval;
          [~,d] = size(U);
          d     = log(d)/log(2);
          Omp0  = 1;
          dQ    = 2;
          
          for indx=1:d
              
             Pplus            = 1/2*(eye(2^(2*d))+ArbDimSWAP(indx,indx+d,2*d));
             Omp0             = Omp0*Pplus;
             
          end
          
          Omp0 = Omp0* 1/( dQ+1 )^d;
          
          ep      = zeros(1,d);
          epNames = cell(1,d);
          
          for indx=1:d
              
              Pminus        = 1/2*(eye(2^(2*d))-ArbDimSWAP(indx,indx+d,2*d));
              ep(indx)      = 2 * trace(  kron(U,U) * Omp0 * kron(U',U') * Pminus) ;
              epNames{indx} = strcat(num2str(indx),'|rest');       
        
              
          end

          obj.one_tangles = {real(ep),epNames};
            

      end
      
      function obj=one_tangles_General_U_Max_Bound(obj)
      %To calculate the one-tangles of general U, based on the partition
      %method of the paper on multipartite entanglement.
      %This bound is not tight! (Only when AME(2n,d) exists).
          U  = obj.Uval;
          D  = log(size(U,2))/log(2);
          dQ = 2;
          
          c  = 1;
          ab = D-1;
        
          expr_bound=0;
        
          for ii=1:D
           
              multip  = nchoosek(D,ii);
              x_prime = ii; 
              y_prime = D-ii;
            
              dim_abx_prime = 2^(ab+x_prime);
              dim_cy_prime  = 2^(c+y_prime);
              min_dim       = min(dim_abx_prime,dim_cy_prime);
              expr_bound    = expr_bound + multip*1/min_dim;
            
          end
          
          x_prime = 0; 
          y_prime = D;
        
          dim_abx_prime = 2^(ab+x_prime);
          dim_cy_prime  = 2^(c+y_prime);
        
          min_dim       = min(dim_abx_prime,dim_cy_prime);
          expr_bound    = expr_bound + 1/min_dim;
          expr_bound    = dQ^D/(dQ+1)^D *expr_bound;
          expr_bound    = 1-expr_bound;
          obj.Max_Bound = expr_bound;
           
      end
      
      function obj=one_tangles_General_U_byparts(obj)
      %To calculate the one-tangles of general U, based on the partition
      %methods of Ref: T. Linowski, G. Rajchel-Mieldzioć, and K. Życzkowski,
      %Entangling Power of Multipartite Unitary Gates, J. Phys. A
      %53, 125303 (2020).
      
          U  = obj.Uval;
          n  = log(size(U,2))/log(2);
          U  = U(:);
          U  = U/norm(U);
          dQ = 2;
           
          byparts_primal    = {1:n};               %Consider only separating 1 subsystem from rest.
          byparts_secondary = all_byparts(n);      %Get all bypartitions of the secondary system
           
          expr   = zeros(1,n);
          labels = cell(1,n);
          cnt    = 0;
          
          for jj=1:length(byparts_primal)  %Loop over primal system bypartitions

              primal_partitions = byparts_primal{jj}; 
              [~,cols_prime]    = size(primal_partitions); 
          
              for kk=1:cols_prime  
              
                  cnt = cnt+1;
                  p   = primal_partitions(:,kk);  

                  %Label the p|q: 
                  labels{cnt} = strcat(num2str(p),'|rest');
              
                  for ii=2:(n+1) %Will add the ii=1 by hand later %loop over secondary system bipartitions
                  
                      xprime       = byparts_secondary{ii};
                      [~,cols_sec] = size(xprime); %dimension of secondary bipartition

                      for ll=1:cols_sec
                      
                          expr(cnt) = expr(cnt) +...
                                      trace(   ptrace( U*U',   [p.' (xprime(:,ll)+n).'] , repmat(2,[1,2*n])  )^2 );  %tracing the p and the x'

                      end
                      
                  end
                  
                  %Add the empty bipartition of the secondary system
                  
                  expr(cnt) = expr(cnt) + trace( ptrace(U*U', p.' , repmat(2,[1,2*n]) )^2  );
                  expr(cnt) = ( 1- dQ^n/(dQ+1)^n * expr(cnt)  );  %Notice here: we drop the 2 factor, following Zanardi's convention.

              end
              
          end
          
          obj.one_tangles={expr,labels};
          
      end
      
  end
  
    
end