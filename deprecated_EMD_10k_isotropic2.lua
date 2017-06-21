arg = { ... }

assert(#arg == 4, "Expected 4 arguments")
assert(argSeed ~= nil, "Expected seed")

prng = DSFMT.create(argSeed)

evolveTime       = arg[1]
reverseOrbitTime = arg[1] / arg[2]

r0  = arg[3]

dwarfMass  = arg[4]

model1Bodies = 100
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.32"

function makePotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
   }
end

encMass = plummerTimestepIntegral(r0, sqr(r0) + sqr(r0/arg[4]) , dwarfMass, 1e-7)

function makeContext()
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = sqr(1/10.0) * sqrt((pi_4_3 * cube(r0)) / (encMass + dwarfMass)),
      eps2       = calculateEps2(totalBodies, r0),
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

    firstModel = predefinedModels.plummer{
        nbody       = model1Bodies,
        prng        = prng,
        position    = finalPosition,
        velocity    = finalVelocity,
        mass        = dwarfMass,
        scaleRadius = r0,
        ignore      = true
    }

      for i,v in ipairs(firstModel)
      do
     	v.ignore = false
      end

return firstModel
end

function makeHistogram()
   return HistogramParams.create()
end


