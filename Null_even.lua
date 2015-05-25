arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)

evolveTime       = arg[1]
reverseOrbitTime = arg[2]

r0  = arg[3]
light_r_ratio = arg[4]

dwarfMass  = arg[5]
light_mass_ratio = arg[6]

model1Bodies = 2000

totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.48"


rscale_l=r0
rscale_d=r0/light_r_ratio
mass_l=dwarfMass * light_mass_ratio
mass_d= dwarfMass*(1.0- light_mass_ratio)


function makePotential()
   
return  nil
end
-- encMass = plummerTimestepIntegral(rscale_l,  rscale_d , mass_d, 1e-7)
-- mass_enc_l= mass_l* (rscale_l)^3* ( (sqr(rscale_l)+ sqr(rscale_l) ) )^(-3.0/2.0)


mass_enc_d= mass_d* (rscale_l)^3* ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)


function makeContext()
   soften_length= (mass_l*rscale_l +mass_d*rscale_d)/(mass_d+mass_l)
   return NBodyCtx.create{
      timeEvolve = evolveTime,
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
    local finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)

    firstModel = predefinedModels.isotropic{
        nbody       = model1Bodies,
        prng        = prng,
        position    = finalPosition,
        velocity    = finalVelocity,
        mass1        = mass_l, --light matter
	mass2       = mass_d,--dark
        scaleRadius1 = rscale_l,
	scaleRadius2 = rscale_d,
        ignore      = true
    }

return firstModel
end

-- function makeHistogram()
--    return HistogramParams.create()
-- end

function makeHistogram()
    return HistogramParams.create{
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = 0.0,
     lambdaEnd = 200,
     lambdaBins = 100,
     betaStart = 0.0,
     betaEnd = -100,
     betaBins = 1
}
end

