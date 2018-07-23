arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)

evolveTime       = tonumber(arg[1])
reverseOrbitTime = arg[1] / tonumber(arg[2])

r0  = tonumber(arg[3])
light_r_ratio = tonumber(arg[4])

dwarfMass  = tonumber(arg[5])
light_mass_ratio = tonumber(arg[6])

model1Bodies = 100000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.42"

function makePotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
   }
end

encMass = plummerTimestepIntegral(r0*light_r_ratio, sqr(r0) + sqr(r0/light_r_ratio) , dwarfMass, 1e-7)

IGNORE = true

if (light_mass_ratio > 0.5) then
   IGNORE = false
end

function makeContext()
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (encMass + dwarfMass)),
      eps2       = calculateEps2(totalBodies, r0/light_r_ratio),
      criterion  = "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end

-- Also required

if (IGNORE) then
   function makeBodies(ctx, potential)
      local firstModel
      local finalPosition, finalVelocity = reverseOrbit{
	 potential = potential,
	 position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
	 velocity  = Vector.create(-156, 79, 107),
	 tstop     = reverseOrbitTime,
	 dt        = ctx.timestep / 10.0
      }
      
      firstModel = predefinedModels.isotropic{
	 nbody       = model1Bodies,
	 prng        = prng,
	 position    = finalPosition,
	 velocity    = finalVelocity,
	 mass1        = dwarfMass * light_mass_ratio,
	 mass2       = dwarfMass - (dwarfMass * light_mass_ratio),
	 scaleRadius1 = r0,
	 scaleRadius2 = r0/light_r_ratio,
	 ignore      = true
      }

      if (tonumber(light_mass_ratio) == 1.0) then
	 for i,v in ipairs(firstModel)
	 do
	    v.ignore = false
	 end
	 
	 
      else 
	 count = 0
	 while (count < light_mass_ratio * totalBodies)
	 do
	    
	    for i,v in ipairs(firstModel)
	    do
	       --Figure out properties
	       r_2 = (finalPosition.x - v.position.x)^2 + (finalPosition.y - v.position.y)^2 + (finalPosition.z - v.position.z)^2
	       r = r_2 ^ (0.5)
	       scale1 = r0 
	       lightDensityToMaxRatio = 1/((1 + r^2/scale1^2)^(5/2))

	       --Get light sphere properties

	       --Chance object is light_mass = light/dwarf
	       chance = lightDensityToMaxRatio
	       
	       ---Do random calculation
	       if(prng:random() > chance and v.ignore==true and count < light_mass_ratio * totalBodies)
	       then
		  v.ignore=false
		  count = count + 1
	       end
	    end
	 end     
      end
  
      return firstModel
   end

else -- if the baryons outnumber the dark matter
   --print("Branch 2\n")
   function makeBodies(ctx, potential)
      local firstModel
      local finalPosition, finalVelocity = reverseOrbit{
         potential = potential,
         position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
         velocity  = Vector.create(-156, 79, 107),
         tstop     = reverseOrbitTime,
         dt        = ctx.timestep / 10.0
      }

      firstModel = predefinedModels.isotropic{
         nbody       = model1Bodies,
         prng        = prng,
         position    = finalPosition,
         velocity    = finalVelocity,
         mass1        = dwarfMass * light_mass_ratio,
         mass2       = dwarfMass - (dwarfMass * light_mass_ratio),
         scaleRadius1 = r0,
         scaleRadius2 = r0/light_r_ratio,
	  ignore      = false
      }

      light_mass_ratio = 1 - light_mass_ratio -- This variable is now poorly named
      if (tonumber(light_mass_ratio) == 1.0) then
         for i,v in ipairs(firstModel)
         do
	        v.ignore = true
         end

      else
         count = 0
         while (count < light_mass_ratio * totalBodies)
         do

            for i,v in ipairs(firstModel)
            do
               --Figure out properties                                                                                                            
               r_2 = (finalPosition.x - v.position.x)^2 + (finalPosition.y - v.position.y)^2 + (finalPosition.z - v.position.z)^2
               r = r_2 ^ (0.5)
               scale1 = r0
	       scale2 = scale1 / light_r_ratio
	       --print(scale2)
               darkDensityToMaxRatio = 1/((1 + r^2/scale2^2)^(5/2))

               --Get dark sphere properties                                                                                                      

               --Chance object is dark mass = light/dwarf
               chance = darkDensityToMaxRatio

               ---Do random calculation                                                                                                           
               if(prng:random() > chance and v.ignore==false and count < light_mass_ratio * totalBodies)
               then
                  v.ignore=true
                  count = count + 1
               end
            end
         end
      end

      return firstModel
   end
end

function makeHistogram()
    return HistogramParams.create{
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = -50,
     lambdaEnd = 50,
     lambdaBins = 50,
     betaStart = -2,
     betaEnd = 2,
     betaBins = 1 
}
end


