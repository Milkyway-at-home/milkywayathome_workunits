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

model1Bodies = 2000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.50"


rscale_l=r0
rscale_d=r0/light_r_ratio
mass_l=dwarfMass * light_mass_ratio
mass_d= dwarfMass*(1.0- light_mass_ratio)


function makePotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
   }
end

-- encMass = plummerTimestepIntegral(r0*light_r_ratio, sqr(r0) + sqr(r0/light_r_ratio) , dwarfMass, 1e-7)
encMass = plummerTimestepIntegral(rscale_l,  rscale_d , mass_d, 1e-7)
mass_enc_d= mass_d* (rscale_l)^3* ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)


function makeContext()
   soften_length= (mass_l*rscale_l +mass_d*rscale_d)/(mass_d+mass_l)
   return NBodyCtx.create{     
      timeEvolve = evolveTime,
--       timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (encMass + dwarfMass)),
      timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(rscale_l)) / (mass_enc_d + mass_l)),
      eps2       = calculateEps2(totalBodies, soften_length),
      criterion  = "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end

-- Also required

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
      mass1       = mass_l,
      mass2       = mass_d,
      scaleRadius1 = rscale_l,
      scaleRadius2 = rscale_d,
      ignore      = true
  }

  return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = -300,
     lambdaEnd = 300,
     lambdaBins = 50,
     betaStart = -20,
     betaEnd = 20,
     betaBins = 1
}
end


