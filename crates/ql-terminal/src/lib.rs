//! A dependency-free terminal dashboard and umbrella SDK for ql-rs.

pub mod prelude {
    pub use crate::TabId;
    pub use ql_ad::*;
    pub use ql_commodities::*;
    pub use ql_commodities::{
        BessAnalyticalSdpDegradationStorageReport, BessAnalyticalSdpDegradationStorageSpec,
        evaluate_bess_analytical_sdp_degradation_storage,
    };
    pub use ql_commodities::{
        WarinBaeuerleStorageRollingHedgeReport, WarinBaeuerleStorageRollingHedgeSpec,
        evaluate_warin_baeuerle_storage_rolling_hedge,
    };
    pub use ql_credit::*;
    pub use ql_credit::{
        BrigoMoriniCmcdsCreditOptionReport, BrigoMoriniCmcdsCreditOptionSpec,
        evaluate_brigo_morini_cmcds_and_credit_index_option,
    };
    pub use ql_credit::{
        DamodaranRsmCountryRiskAndDebtModReport, DamodaranRsmCountryRiskAndDebtModSpec,
        FsbFinmaGsibResolutionBailInWaterfallReport, FsbFinmaGsibResolutionBailInWaterfallSpec,
        evaluate_damodaran_rsm_country_risk_and_debt_modification,
        evaluate_fsb_finma_gsib_resolution_bail_in_waterfall,
    };
    pub use ql_credit::{
        FedsDeSpiegeleerAt1CocoDualTriggerReport, FedsDeSpiegeleerAt1CocoDualTriggerSpec,
        evaluate_feds_despiegeleer_at1_coco_dual_trigger,
    };
    pub use ql_credit::{
        LandoIsraelGeneratorLogarithmReport, LandoIsraelGeneratorLogarithmSpec,
        compute_lando_israel_generator_logarithm_and_migration,
    };
    pub use ql_credit::{
        RsmAsc470DebtModificationReport, RsmAsc470DebtModificationSpec,
        evaluate_rsm_asc470_debt_modification_and_restructuring,
    };
    pub use ql_curves::*;
    pub use ql_curves::{
        AdbRbiNdirsAssetSwapReport, AdbRbiNdirsAssetSwapSpec,
        evaluate_adb_rbi_ndirs_and_asset_swap_basis,
    };
    pub use ql_curves::{
        BcbsD368SixScenarioIrrbbReport, BcbsD368SixScenarioIrrbbSpec,
        evaluate_bcbs_d368_six_scenario_irrbb,
    };
    pub use ql_curves::{
        FabozziKeyRateLscImmunizationReport, FabozziKeyRateLscImmunizationSpec,
        evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization,
    };
    pub use ql_curves::{
        RomerVeronesiTaylorRuleMacroCurveReport, RomerVeronesiTaylorRuleMacroCurveSpec,
        evaluate_romer_veronesi_taylor_rule_macro_yield_curve,
    };
    pub use ql_dml::*;
    pub use ql_dml::{
        AgrawalBoydDifferentiableQpReport, AgrawalBoydDifferentiableQpSpec,
        solve_agrawal_boyd_differentiable_qp_layer,
    };
    pub use ql_dml::{
        GlassermanSaddlepointIsStratifiedLhsReport, GlassermanSaddlepointIsStratifiedLhsSpec,
        evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad,
    };
    pub use ql_dml::{
        GreeneWooldridgeHeckmanSelectionGmmReport, GreeneWooldridgeHeckmanSelectionGmmSpec,
        HastieTibshiraniElasticNetGcvReport, HastieTibshiraniElasticNetGcvSpec,
        evaluate_greene_wooldridge_heckman_survivorship_and_gmm,
        solve_hastie_tibshirani_elastic_net_and_gcv,
    };
    pub use ql_dml::{
        HugeSavineDifferentialPcaReport, HugeSavineDifferentialPcaSpec,
        evaluate_huge_savine_differential_pca_and_twin_dml,
    };
    pub use ql_exotics::*;
    pub use ql_exotics::{
        ChoiLiDengMultiAssetSpreadBasketReport, ChoiLiDengMultiAssetSpreadBasketSpec,
        evaluate_choi_li_deng_multi_asset_spread_basket,
    };
    pub use ql_exotics::{
        CurranRogersShiArithmeticAsianReport, CurranRogersShiArithmeticAsianSpec,
        price_curran_rogers_shi_arithmetic_asian,
    };
    pub use ql_exotics::{
        SheppShiryaevShrevePerpetualRussianReport, SheppShiryaevShrevePerpetualRussianSpec,
        price_shepp_shiryaev_shreve_perpetual_russian_and_mckean,
    };
    pub use ql_fixed_income::*;
    pub use ql_fixed_income::{
        AcharyaFedsBankRunAndCocoReport, AcharyaFedsBankRunAndCocoSpec,
        analyze_acharya_feds_bank_run_and_optimal_coco,
    };
    pub use ql_fixed_income::{
        AnDaoquanCliffordChanceResettableCbReport, AnDaoquanCliffordChanceResettableCbSpec,
        evaluate_an_daoquan_clifford_chance_resettable_convertible,
    };
    pub use ql_fixed_income::{
        AsquithSirbuShreveCbGameReport, AsquithSirbuShreveCbGameSpec,
        DavidsonLevinCmoPacCompanionReport, DavidsonLevinCmoPacCompanionSpec,
        FasbAsu2024InducedConversionCpdiReport, FasbAsu2024InducedConversionCpdiSpec,
        evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game,
        evaluate_davidson_levin_cmo_pac_companion_waterfall,
        evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid,
    };
    pub use ql_fixed_income::{
        KoziolAsquithFasbStrategicConvertibleReport, KoziolAsquithFasbStrategicConvertibleSpec,
        VeronesiFabozziConvexityPrepaymentReport, VeronesiFabozziConvexityPrepaymentSpec,
        evaluate_koziol_asquith_fasb_strategic_convertible,
        evaluate_veronesi_fabozzi_futures_convexity_and_prepayment,
    };
    pub use ql_fixed_income::{
        MayerBrownSsgaCappedCallMakeWholeReport, MayerBrownSsgaCappedCallMakeWholeSpec,
        evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole,
    };
    pub use ql_fx::*;
    pub use ql_fx::{
        BaranBianchettiMtmXccyBasisReport, BaranBianchettiMtmXccyBasisSpec,
        evaluate_baran_bianchetti_mtm_xccy_basis_swap,
    };
    pub use ql_fx::{
        CrossCurrencyConvertibleWarrantReport, CrossCurrencyConvertibleWarrantSpec,
        price_cross_currency_convertible_and_warrant_bond,
    };
    pub use ql_inflation::*;
    pub use ql_inflation::{
        CochraneFiscalTheoryPriceLevelReport, CochraneFiscalTheoryPriceLevelSpec,
        evaluate_cochrane_ftpl_and_maturity_inflation_smoothing,
    };
    pub use ql_instruments::*;
    pub use ql_levy::{
        CgmyFractionalFftAndPideReport, CgmyFractionalFftAndPideSpec, CgmyProcess,
        ContTankovBatesSvjPideReport, ContTankovBatesSvjPideSpec, KouJumpDiffusion, LevyModel,
        MertonJumpDiffusion, VarianceGamma, cos_levy_price, cos_levy_price_strip,
        evaluate_cgmy_fractional_fft_and_pide,
        evaluate_cont_tankov_bates_svj_characteristic_and_pide, merton_series_price,
    };
    pub use ql_levy::{
        ContTankovMertonJumpDiffusionReport, ContTankovMertonJumpDiffusionSpec,
        carr_madan_fft_levy_price_strip, evaluate_cont_tankov_merton_jump_diffusion,
        lewis_levy_price_strip,
    };
    pub use ql_math::*;
    pub use ql_mc::*;
    pub use ql_models::{
        AcmCochranePiazzesiTermPremiumReport, AcmCochranePiazzesiTermPremiumSpec,
        BansalYaronLongRunRiskReport, BansalYaronLongRunRiskSpec, BansalYaronLongRunRisksReport,
        BansalYaronLongRunRisksSpec, BernankeGertlerGilchristAcceleratorReport,
        BernankeGertlerGilchristAcceleratorSpec, BlanchardLhuillierLorenzoniNewsNoiseReport,
        BlanchardLhuillierLorenzoniNewsNoiseSpec, BroadieGlassermanKouDiscreteBarrierReport,
        BroadieGlassermanKouDiscreteBarrierSpec, BrunnermeierSannikovVolatilityParadoxAnalytics,
        BrunnermeierSannikovVolatilityParadoxSpec, CampbellCochraneExternalHabitReport,
        CampbellCochraneExternalHabitSpec, CampbellShillerCapeExcessYieldReport,
        CampbellShillerCapeExcessYieldSpec, CampbellVuolteenahoBadGoodBetaIcapmAnalytics,
        CampbellVuolteenahoBadGoodBetaIcapmSpec, ChristianoEichenbaumEvansNkAnalytics,
        ChristianoEichenbaumEvansNkSpec, DiamondDybvigGoldsteinPauznerBankRunAnalytics,
        DiamondDybvigGoldsteinPauznerBankRunSpec, DrechslerSavovSchnablDepositChannelAnalytics,
        DrechslerSavovSchnablDepositChannelSpec, EpsteinZinWeilRecursiveUtilitySdfReport,
        EpsteinZinWeilRecursiveUtilitySdfSpec, FarhiGabaixRareDisasterExchangeRatePuzzleAnalytics,
        FarhiGabaixRareDisasterExchangeRatePuzzleSpec, GabaixRareDisasterResilienceReport,
        GabaixRareDisasterResilienceSpec, GabaixRareDisastersResilienceAnalytics,
        GabaixRareDisastersResilienceSpec, GabaixVariableRareDisasterReport,
        GabaixVariableRareDisasterSpec, GertlerKaradiUnconventionalMonetaryPolicyAnalytics,
        GertlerKaradiUnconventionalMonetaryPolicySpec, HeKrishnamurthyIntermediarySdfAnalytics,
        HeKrishnamurthyIntermediarySdfSpec, JermannQuadriniFinancialShockAnalytics,
        JermannQuadriniFinancialShockSpec, KiyotakiMooreCollateralCycleReport,
        KiyotakiMooreCollateralCycleSpec, KollerMckinseyKvdEvaValuationReport,
        KollerMckinseyKvdEvaValuationSpec, LaubachWilliamsNaturalRateRStarReport,
        LaubachWilliamsNaturalRateRStarSpec, LettauLudvigsonCayConditionalRiskPremiumAnalytics,
        LettauLudvigsonCayConditionalRiskPremiumSpec, LettauLudvigsonCayForecastReport,
        LettauLudvigsonCayForecastSpec, LjungqvistNewKeynesianReeReport,
        LjungqvistNewKeynesianReeSpec, LjungqvistSargentMccallSearchReport,
        LjungqvistSargentMccallSearchSpec, LucasMehraPrescottEndowmentTreeReport,
        LucasMehraPrescottEndowmentTreeSpec, OhlsonPenmanLinearInfoDynamicsReport,
        OhlsonPenmanLinearInfoDynamicsSpec, PastorVeronesiPoliticalUncertaintyAnalytics,
        PastorVeronesiPoliticalUncertaintySpec, PenmanNissimDuPontReoiValuationReport,
        PenmanNissimDuPontReoiValuationSpec, PenmanReoiOhlsonJuettnerReport,
        PenmanReoiOhlsonJuettnerSpec, PiazzesiSchneiderNominalBondAnalytics,
        PiazzesiSchneiderNominalBondSpec, PiazzesiSwansonFomcFuturesHikeReport,
        PiazzesiSwansonFomcFuturesHikeSpec, WachterDisasterVolatilityAnalytics,
        WachterDisasterVolatilitySpec, evaluate_acm_cochrane_piazzesi_affine_term_premium,
        evaluate_bansal_yaron_long_run_risk_stochastic_vol_equity_premium,
        evaluate_bansal_yaron_long_run_risks_sdf,
        evaluate_bernanke_gertler_gilchrist_financial_accelerator,
        evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter,
        evaluate_broadie_glasserman_kou_discrete_barrier_lookback,
        evaluate_brunnermeier_sannikov_intermediary_volatility_paradox,
        evaluate_campbell_cochrane_external_habit_sdf,
        evaluate_campbell_shiller_cape_excess_yield_forecast,
        evaluate_campbell_vuolteenaho_bad_beta_good_beta_icapm,
        evaluate_christiano_eichenbaum_evans_calvo_investment_nk_impulse,
        evaluate_diamond_dybvig_goldstein_pauzner_bank_run_global_game,
        evaluate_drechsler_savov_schnabl_deposit_channel_liquidity_premium,
        evaluate_epstein_zin_weil_recursive_utility_sdf,
        evaluate_farhi_gabaix_rare_disaster_exchange_rate_uip,
        evaluate_gabaix_rare_disaster_resilience,
        evaluate_gabaix_rare_disasters_time_varying_resilience_equity_premium,
        evaluate_gabaix_variable_rare_disaster_resilience_equity_puzzle,
        evaluate_gertler_karadi_intermediary_credit_policy,
        evaluate_he_krishnamurthy_intermediary_capital_constraint_sdf,
        evaluate_jermann_quadrini_financial_shocks_enforcement_constraint_credit_cycle,
        evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier,
        evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation,
        evaluate_laubach_williams_natural_rate_r_star,
        evaluate_lettau_ludvigson_cay_conditional_equity_premium,
        evaluate_lettau_ludvigson_cay_consumption_wealth_forecast,
        evaluate_ljungqvist_new_keynesian_three_equation_ree,
        evaluate_ljungqvist_sargent_mccall_bellman_search,
        evaluate_lucas_mehra_prescott_endowment_tree_equity_premium,
        evaluate_ohlson_penman_linear_information_dynamics,
        evaluate_pastor_veronesi_political_uncertainty_risk_premium,
        evaluate_penman_nissim_dupont_and_reoi_valuation,
        evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation,
        evaluate_piazzesi_schneider_stagflation_nominal_bond_term_premium,
        evaluate_piazzesi_swanson_fomc_futures_hike_tree,
        evaluate_wachter_time_varying_disaster_intensity_equity_volatility,
    };
    pub use ql_models::{
        AnalyticBsmHullWhiteEngine, BatesParams, Bergomi2FParams,
        BergomiForwardStartAndVarStripReport, BergomiForwardStartAndVarStripSpec,
        CalibratedHestonModel, Complex, CosHestonEngine, CosHestonHullWhiteEngine,
        DamodaranRdCapitalizationRealOptionReport, DamodaranRdCapitalizationRealOptionSpec,
        DoubleHestonParams, HestonHullWhiteParams, HestonMarketQuote, HestonParams, HullWhite1F,
        LiftedHestonParams, MasColellCochraneArrowDebreuSdfReport,
        MasColellCochraneArrowDebreuSdfSpec, MixedRoughBergomiVixReport, MixedRoughBergomiVixSpec,
        QuinticOuParams, QuinticOuVixSpxReport, VixOptionValuation, bergomi_2f_implied_vol_bg,
        bergomi_2f_mc_price, bergomi_2f_smile_moments_bg, calibrate_heston,
        double_heston_cos_price, double_heston_cos_strip_prices,
        evaluate_bergomi_forward_start_and_var_strip,
        evaluate_damodaran_rd_capitalization_and_patent_real_option,
        evaluate_mascolell_cochrane_arrow_debreu_sdf, evaluate_mixed_rough_bergomi_vix,
        evaluate_quintic_ou_spx_vix, h1_hw_char_func_unit_spot,
        heston_bates_vix_future_and_option_price, heston_cos_price, lifted_heston_char_fn,
        lifted_heston_char_func_unit_spot, lifted_heston_cos_price, lifted_heston_kernel,
        local_vol::{
            AndreasenHugeGridSpec, AndreasenHugeSliceResult, FeldmanSabrCalibrationResult,
            ZabrParams, ZabrSmileQuote, evaluate_zabr_smile, feldman_analytic_sabr_calibration,
            interpolate_andreasen_huge_call, solve_andreasen_huge_step,
        },
        replicate_cboe_vix_index,
    };
    pub use ql_models::{
        CarmonaDelarueMeanFieldSystemicReport, CarmonaDelarueMeanFieldSystemicSpec,
        solve_carmona_delarue_mean_field_systemic_interbank_game,
    };
    pub use ql_models::{
        DoubleHestonKernelSlvReport, DoubleHestonKernelSlvSpec,
        evaluate_double_heston_cumulants_and_kernel_slv,
    };
    pub use ql_models::{
        GairatPolalaCheyetteCapletReport, GairatPolalaCheyetteCapletSpec,
        GatheralJusselinQuadraticRoughHestonReport, GatheralJusselinQuadraticRoughHestonSpec,
        SchilitSloanDechowDichevQualityReport, SchilitSloanDechowDichevQualitySpec,
        evaluate_gairat_polala_cheyette_local_vol_caplet,
        evaluate_gatheral_jusselin_quadratic_rough_heston,
        evaluate_schilit_sloan_dechow_dichev_earnings_quality,
    };
    pub use ql_models::{
        JacquierDeepRoughVolSurrogateReport, JacquierDeepRoughVolSurrogateSpec,
        evaluate_jacquier_chebyshev_rough_vol_surrogate,
    };
    pub use ql_models::{rough::*, slv::*};
    pub use ql_pde::*;
    pub use ql_pde::{
        HaentjensInthoutHhwAdiReport, HaentjensInthoutHhwAdiSpec,
        evaluate_haentjens_inthout_heston_hull_white_adi,
    };
    pub use ql_portfolio::*;
    pub use ql_portfolio::{
        AhmadiJavidEntropicVarPortfolioReport, AhmadiJavidEntropicVarPortfolioSpec,
        CottonSchurHrpRmtNcoReport, CottonSchurHrpRmtNcoSpec, TsayBekkGarchSpilloverReport,
        TsayBekkGarchSpilloverSpec, evaluate_ahmadi_javid_entropic_var_portfolio,
        evaluate_cotton_schur_hrp_and_rmt_nco, evaluate_tsay_bekk_garch_volatility_spillover,
    };
    pub use ql_portfolio::{
        BarberisJinWangProspectAnomalyReport, BarberisJinWangProspectAnomalySpec,
        ChristoffersenBouchaudVarBacktestAndCumulantReport,
        ChristoffersenBouchaudVarBacktestAndCumulantSpec,
        ColinBaconCampisiFixedIncomeAttributionReport, ColinBaconCampisiFixedIncomeAttributionSpec,
        FoucaultVanDerMerweLobLiquidityLvarReport, FoucaultVanDerMerweLobLiquidityLvarSpec,
        HartMasColellRegretMatchingPortfolioReport, HartMasColellRegretMatchingPortfolioSpec,
        evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly,
        evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es,
        evaluate_colin_bacon_campisi_fixed_income_attribution,
        evaluate_foucault_vandermerwe_lob_schedule_and_lvar,
        evaluate_hart_mascolell_regret_matching_allocation,
    };
    pub use ql_portfolio::{
        BouchaudPottersCumulantSmileReport, BouchaudPottersCumulantSmileSpec,
        BrockHommesAdaptiveBeliefReport, BrockHommesAdaptiveBeliefSpec,
        DgtwTreynorMazuyTimingReport, DgtwTreynorMazuyTimingSpec,
        EmbrechtsHillPickandsSubexponentialReport, EmbrechtsHillPickandsSubexponentialSpec,
        KarnoskySingerMencheroAttributionReport, KarnoskySingerMencheroAttributionSpec,
        evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling,
        evaluate_brock_hommes_adaptive_belief_switching,
        evaluate_dgtw_characteristic_and_treynor_mazuy_timing,
        evaluate_embrechts_hill_pickands_and_subexponential_tail,
        evaluate_karnosky_singer_multicurrency_and_menchero_attribution,
    };
    pub use ql_portfolio::{
        BussetiBoydRiskConstrainedKellyReport, BussetiBoydRiskConstrainedKellySpec,
        HarrisRollCorwinSchultzMicrostructureReport, HarrisRollCorwinSchultzMicrostructureSpec,
        analyze_harris_roll_corwin_schultz_microstructure,
        solve_busseti_boyd_risk_constrained_kelly,
    };
    pub use ql_portfolio::{
        CarteaJaimungalPenalvaLobFillHjbReport, CarteaJaimungalPenalvaLobFillHjbSpec,
        GrinoldQianChristophersonAlphaDecayReport, GrinoldQianChristophersonAlphaDecaySpec,
        MadhavanAbnerHillEtfApArbitrageReport, MadhavanAbnerHillEtfApArbitrageSpec,
        evaluate_grinold_qian_christopherson_ic_decay_and_timing,
        evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha,
        solve_cartea_jaimungal_penalva_lob_fill_hjb,
    };
    pub use ql_portfolio::{
        FollmerSchiedGennaioliDiagnosticRiskReport, FollmerSchiedGennaioliDiagnosticRiskSpec,
        MitchellPedersenKratzDarkPoolSpiralReport, MitchellPedersenKratzDarkPoolSpiralSpec,
        evaluate_follmer_schied_convex_risk_and_diagnostic_expectations,
        evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool,
    };
    pub use ql_rates::*;
    pub use ql_rates::{
        BrigoMercurioG2TwistSwaptionReport, BrigoMercurioG2TwistSwaptionSpec,
        price_brigo_mercurio_g2_twist_and_swaption,
    };
    pub use ql_rates::{
        TurfusWillemsRfrInArrearsSmileReport, TurfusWillemsRfrInArrearsSmileSpec,
        evaluate_turfus_willems_rfr_in_arrears_caplet_smile,
    };
    pub use ql_scripting::*;
    pub use ql_time::*;
    pub use ql_vol::*;
    pub use ql_vol::{
        AndreasenHugeOneStepZabrReport, AndreasenHugeOneStepZabrSpec,
        evaluate_andreasen_huge_one_step_zabr_calibration,
    };
    pub use ql_vol::{
        BourgeyGatheralBergomiGuyonSmileReport, BourgeyGatheralBergomiGuyonSmileSpec,
        evaluate_bourgey_gatheral_bergomi_guyon_expansion,
    };
    pub use ql_vol::{
        GatheralDemeterfiCorridorVarianceSkewReport, GatheralDemeterfiCorridorVarianceSkewSpec,
        evaluate_gatheral_demeterfi_corridor_variance_and_skew,
    };
    pub use ql_vol::{
        MingoneGlobalArbitrageFreeEssviReport, MingoneGlobalArbitrageFreeEssviSpec,
        evaluate_mingone_global_arbitrage_free_essvi,
    };
    pub use ql_xva::*;
    pub use ql_xva::{
        AlbaneseCrepeyMporBalanceSheetReport, AlbaneseCrepeyMporBalanceSheetSpec,
        compute_albanese_crepey_mpor_and_balance_sheet_xva,
    };
    pub use ql_xva::{
        AlbaneseCrepeyZeronBalanceSheetXvaReport, AlbaneseCrepeyZeronBalanceSheetXvaSpec,
        evaluate_albanese_crepey_zeron_balance_sheet_xva,
    };
    pub use ql_xva::{
        BcbsD424OutputFloorAndSmaOpRiskReport, BcbsD424OutputFloorAndSmaOpRiskSpec,
        evaluate_bcbs_d424_output_floor_and_sma_op_risk,
    };
    pub use ql_xva::{
        GregoryCsaThresholdMtaCtdColvaReport, GregoryCsaThresholdMtaCtdColvaSpec,
        evaluate_gregory_csa_threshold_mta_and_ctd_colva,
    };
    pub use ql_xva::{
        TurlakovGreenKenyonWwrKvaReport, TurlakovGreenKenyonWwrKvaSpec,
        evaluate_turlakov_green_kenyon_wwr_and_kva,
    };
}

use ql_ad::Real;
use ql_commodities::{
    BessAnalyticalSdpDegradationStorageSpec, evaluate_bess_analytical_sdp_degradation_storage,
};
use ql_commodities::{
    BessSpec, CleanSparkDarkSpreadSpec, ClewlowStricklandSwaptionSpec, CommodityFuturesCurve,
    DegreeDayKind, GasStorageSpec, SchwartzSmith2F, SwingOptionSpec, TakeOrPayContractSpec,
    WeatherDerivativeSpec, kirk_spread_option, price_bess_arbitrage_dp,
    price_clean_spark_dark_spread, price_clewlow_strickland_strip_swaption, price_gas_storage,
    price_swing_option, price_take_or_pay_contract, price_weather_hdd_cdd_option,
};
#[allow(unused_imports)]
use ql_commodities::{
    EydelandWolyniecRampOutageVtaReport, EydelandWolyniecRampOutageVtaSpec,
    evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling,
};
#[allow(unused_imports)]
use ql_commodities::{
    WarinBaeuerleStorageRollingHedgeReport, WarinBaeuerleStorageRollingHedgeSpec,
    evaluate_warin_baeuerle_storage_rolling_hedge,
};
use ql_credit::{
    BaseCorrelationTrancheQuote, BlackCoxStructuralSpec, CalibratedCreditCurve, CdsQuote,
    CmCdsMarketModelSpec, CreditDefaultSwap, IsdaCdsEngine, JltRatingMigrationSpec,
    OKaneBespokeBaseCorrelationSpec, ProtectionSide, SurvivalProbabilityCurve,
    analyze_structural_credit_black_cox, bootstrap_base_correlation_curve,
    evaluate_jlt_rating_migration_credit, price_cmcds_market_model,
    price_okane_bespoke_tranche_and_student_t_basket,
};
use ql_credit::{
    FedsDeSpiegeleerAt1CocoDualTriggerSpec, evaluate_feds_despiegeleer_at1_coco_dual_trigger,
};
#[allow(unused_imports)]
use ql_credit::{
    LandoIsraelGeneratorLogarithmReport, LandoIsraelGeneratorLogarithmSpec,
    compute_lando_israel_generator_logarithm_and_migration,
};
#[allow(unused_imports)]
use ql_credit::{
    SchonbucherDefaultableRangeAccrualCapReport, SchonbucherDefaultableRangeAccrualCapSpec,
    evaluate_schonbucher_defaultable_range_accrual_and_spread_cap,
};
use ql_curves::{BcbsD368SixScenarioIrrbbSpec, evaluate_bcbs_d368_six_scenario_irrbb};
use ql_curves::{DiscountCurve, FlatForwardCurve, MonotoneConvexZeroCurve, YieldTermStructure};
#[allow(unused_imports)]
use ql_dml::{
    AgrawalBoydDifferentiableQpReport, AgrawalBoydDifferentiableQpSpec,
    solve_agrawal_boyd_differentiable_qp_layer,
};
#[allow(unused_imports)]
use ql_dml::{
    GilesCapriottiMultilevelMcAadReport, GilesCapriottiMultilevelMcAadSpec,
    evaluate_giles_capriotti_multilevel_mc_aad,
};
#[allow(unused_imports)]
use ql_dml::{
    GlassermanSaddlepointIsStratifiedLhsReport, GlassermanSaddlepointIsStratifiedLhsSpec,
    evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad,
};
use ql_dml::{HugeSavineDifferentialPcaSpec, evaluate_huge_savine_differential_pca_and_twin_dml};
use ql_exotics::{
    AtlasOptionSpec, BarrierOptionSpec, BarrierType as ExoticBarrierType, EverestOptionSpec,
    ExcursionClockKind, FiniteCapTimerOptionSpec, GammaSwapAndVarOptionSpec, MonitoringFrequency,
    ParisianBarrierOptionSpec, VarianceSwapStrikeQuote, WorstOfAutocallableSpec,
    price_atlas_option_mc_ad, price_barrier_option, price_capped_timer_option_quadrature,
    price_everest_option_mc_ad, price_gamma_swap_and_variance_option, price_parisian_option_mc_ad,
    price_worst_of_autocallable,
};
use ql_exotics::{
    ChoiLiDengMultiAssetSpreadBasketSpec, evaluate_choi_li_deng_multi_asset_spread_basket,
};
#[allow(unused_imports)]
use ql_exotics::{
    CurranRogersShiArithmeticAsianReport, CurranRogersShiArithmeticAsianSpec,
    price_curran_rogers_shi_arithmetic_asian,
};
#[allow(unused_imports)]
use ql_exotics::{
    HeynenKatPartialFractionalLookbackReport, HeynenKatPartialFractionalLookbackSpec,
    evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread,
};
#[allow(unused_imports)]
use ql_exotics::{
    SheppShiryaevShrevePerpetualRussianReport, SheppShiryaevShrevePerpetualRussianSpec,
    price_shepp_shiryaev_shreve_perpetual_russian_and_mckean,
};
#[allow(unused_imports)]
use ql_fixed_income::{
    AcharyaFedsBankRunAndCocoReport, AcharyaFedsBankRunAndCocoSpec,
    analyze_acharya_feds_bank_run_and_optimal_coco,
};
use ql_fixed_income::{
    AfvConvertibleBondSpec, BondAssetSwapSpec, CallableBondOasSpec, CmoTrancheSpec,
    CmoWaterfallSpec, CocoBondSpec, CocoTriggerPayoff, ConvertibleBondSpec,
    DeSpiegeleerCoCoMandatoryConvSpec, MbsBurnoutOasSpec, NelsonSiegelSvensson,
    SpecifiedPoolMbsSpec, TbaDollarRollSpec, TreasuryFutureDeliverable, TreasuryFuturesBasketSpec,
    TuckmanKeyRateDurationSpec, analyze_bond_asset_swap_and_z_spread, analyze_callable_bond_oas_ad,
    analyze_mbs_burnout_oas_ad, analyze_tba_dollar_roll, analyze_treasury_futures_ctd,
    analyze_tuckman_key_rate_durations_and_pca_hedge, price_cmo_waterfall, price_coco_bond_ds,
    price_convertible_bond_afv_pde, price_convertible_bond_tf,
    price_despiegeleer_coco_and_mandatory_convertible, price_specified_pool_mbs,
};
#[allow(unused_imports)]
use ql_fixed_income::{
    KoziolAsquithFasbStrategicConvertibleReport, KoziolAsquithFasbStrategicConvertibleSpec,
    VeronesiFabozziConvexityPrepaymentReport, VeronesiFabozziConvexityPrepaymentSpec,
    evaluate_koziol_asquith_fasb_strategic_convertible,
    evaluate_veronesi_fabozzi_futures_convexity_and_prepayment,
};
#[allow(unused_imports)]
use ql_fixed_income::{
    TuckmanBurghardtTreasuryWildcardBasisReport, TuckmanBurghardtTreasuryWildcardBasisSpec,
    evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch,
};
use ql_fx::{BaranBianchettiMtmXccyBasisSpec, evaluate_baran_bianchetti_mtm_xccy_basis_swap};
use ql_fx::{
    ConstNotionalXccs, FxAccumulatorSpec, FxKikoOptionSpec, FxSmileQuote, FxTarfCapConvention,
    FxTarfSpec, FxVannaVolgaSmile, FxWindowBarrierSpec, MtmResettingXccs, PrdcNoteSpec,
    WystupFxStructuredProductSpec, evaluate_wystup_fx_accumulator_and_dcd,
    garman_kohlhagen_spot_delta, price_fx_accumulator_strip_no_ko, price_fx_kiko_option,
    price_fx_tarf_mc_ad, price_fx_window_barrier_option, price_mtm_resetting_xccs,
    price_prdc_note_hybrid,
};
#[allow(unused_imports)]
use ql_fx::{
    CrossCurrencyConvertibleWarrantReport, CrossCurrencyConvertibleWarrantSpec,
    price_cross_currency_convertible_and_warrant_bond,
};
#[allow(unused_imports)]
use ql_fx::{
    WystupPivotTarfCorridorReport, WystupPivotTarfCorridorSpec,
    evaluate_wystup_pivot_tarf_and_corridor_strip,
};
use ql_inflation::{
    BrigoMercurioInflationMarketModelSpec, CalibratedInflationCurve, JarrowYildirimParams,
    MultiplicativeSeasonality, YoYInflationCapFloor, YoYInflationCapFloorType, ZciisCapFloorSpec,
    ZeroCouponInflationSwap, ZeroInflationCurve, price_brigo_mercurio_yoy_and_lpi_market_model,
    price_yoy_capfloor_jarrow_yildirim, price_zciis_capfloor_jarrow_yildirim,
};
#[allow(unused_imports)]
use ql_inflation::{
    CochraneFiscalTheoryPriceLevelReport, CochraneFiscalTheoryPriceLevelSpec,
    evaluate_cochrane_ftpl_and_maturity_inflation_smoothing,
};
use ql_instruments::{BaroneAdesiWhaleyEngine, OptionType, SwapType, VanillaSwapSpec};
use ql_levy::{CgmyProcess, VarianceGamma, cos_levy_price, evaluate_cgmy_fractional_fft_and_pide};
#[allow(unused_imports)]
use ql_levy::{
    ContTankovBatesSvjPideReport, ContTankovBatesSvjPideSpec,
    evaluate_cont_tankov_bates_svj_characteristic_and_pide,
};
use ql_models::rough::{
    QuadraticRoughHestonParams, RoughBergomiParams, qrh_implied_vol_expansion,
    rbergomi_hybrid_mc_price,
};
use ql_models::{
    AndreasenHugeGridSpec, BatesParams, Bergomi2FParams, DoubleHestonParams, HestonParams,
    HullWhite1F, LiftedHestonParams, MixedRoughBergomiVixSpec, QuinticOuParams,
    bergomi_2f_implied_vol_bg, double_heston_cos_price, evaluate_mixed_rough_bergomi_vix,
    evaluate_quintic_ou_spx_vix, heston_bates_vix_future_and_option_price,
    interpolate_andreasen_huge_call, lifted_heston_cos_price,
    local_vol::{ZabrParams, evaluate_zabr_smile},
    rough::{GuyonPdvParams, pdv_implied_vol_expansion},
    solve_andreasen_huge_step,
};
#[allow(unused_imports)]
use ql_models::{
    CarmonaDelarueMeanFieldSystemicReport, CarmonaDelarueMeanFieldSystemicSpec,
    solve_carmona_delarue_mean_field_systemic_interbank_game,
};
#[allow(unused_imports)]
use ql_models::{
    DoubleHestonKernelSlvReport, DoubleHestonKernelSlvSpec,
    evaluate_double_heston_cumulants_and_kernel_slv,
};
use ql_models::{
    JacquierDeepRoughVolSurrogateSpec, evaluate_jacquier_chebyshev_rough_vol_surrogate,
};
use ql_portfolio::{
    AlfonsiFruthSchiedLobSpec, BoxConstrainedMarkowitzSpec, DarkLitRoutingSpec, EntropyPoolingSpec,
    ExpectileRiskAndPortfolioSpec, IdzorekBlackLittermanSpec, KellyAllocationSpec, MarketSnapshot,
    ObizhaevaWangSpec, PortfolioDag, PortfolioRiskSummary, PortfolioTrade, RobustMarkowitzSpec,
    box_constrained_markowitz_weights, effective_number_of_bets_minimum_torsion,
    entropy_pooling_posterior, equal_risk_contribution_weights,
    evaluate_and_optimize_expectile_portfolio, gerber_covariance_matrix,
    hierarchical_risk_parity_weights, idzorek_black_litterman_posterior,
    kelly_fractional_allocation, nested_clustered_optimization_weights, obizhaeva_wang_execution,
    optimal_dark_lit_execution, robust_markowitz_weights, solve_alfonsi_fruth_schied_lob_execution,
};
use ql_portfolio::{
    AvellanedaStoikovMmSpec, MeanCdarPortfolioSpec, SchurHrpSpec, conditional_drawdown_at_risk,
    entropic_value_at_risk, factor_risk_attribution, gaussian_expected_shortfall,
    oas_and_nonlinear_shrinkage_covariance, optimize_mean_cdar_portfolio, schur_complementary_hrp,
    solve_avellaneda_stoikov_market_making,
};
#[allow(unused_imports)]
use ql_portfolio::{
    BussetiBoydRiskConstrainedKellyReport, BussetiBoydRiskConstrainedKellySpec,
    HarrisRollCorwinSchultzMicrostructureReport, HarrisRollCorwinSchultzMicrostructureSpec,
    analyze_harris_roll_corwin_schultz_microstructure, solve_busseti_boyd_risk_constrained_kelly,
};
#[allow(unused_imports)]
use ql_portfolio::{
    CarteaJaimungalPenalvaLobFillHjbReport, CarteaJaimungalPenalvaLobFillHjbSpec,
    GrinoldQianChristophersonAlphaDecayReport, GrinoldQianChristophersonAlphaDecaySpec,
    MadhavanAbnerHillEtfApArbitrageReport, MadhavanAbnerHillEtfApArbitrageSpec,
    evaluate_grinold_qian_christopherson_ic_decay_and_timing,
    evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha,
    solve_cartea_jaimungal_penalva_lob_fill_hjb,
};
#[allow(unused_imports)]
use ql_portfolio::{
    EngleRussellAcdTradeDurationReport, EngleRussellAcdTradeDurationSpec,
    GrinoldKahnClarkeFundamentalLawReport, GrinoldKahnClarkeFundamentalLawSpec,
    HarveyLiuZhuMultipleTestingSharpeReport, HarveyLiuZhuMultipleTestingSharpeSpec,
    evaluate_engle_russell_acd_trade_duration_and_hazard,
    evaluate_grinold_kahn_clarke_fundamental_law,
    evaluate_harvey_liu_zhu_multiple_testing_and_mintrl,
};
use ql_portfolio::{
    FollmerSchiedGennaioliDiagnosticRiskSpec, MitchellPedersenKratzDarkPoolSpiralSpec,
    evaluate_follmer_schied_convex_risk_and_diagnostic_expectations,
    evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool,
};
use ql_rates::{
    BermudanSwaption, CheyetteSwaptionSpec, EuropeanSwaption, G2PlusPlus,
    HullWhiteBermudanTreeEngine, MarkovFunctionalSwaptionSpec, RfrCapletSpec, RfrRateType,
    SinhOuRfrCapletSpec, SwaptionSettlement, TenorBasisSwapSpec,
    price_cheyette_swaption_markov_projection, price_markov_functional_1f,
    price_sinh_ou_rfr_caplet, price_tenor_basis_swap, willems_backward_sabr_caplet_price,
};
#[allow(unused_imports)]
use ql_rates::{
    BrigoMercurioG2TwistSwaptionReport, BrigoMercurioG2TwistSwaptionSpec,
    price_brigo_mercurio_g2_twist_and_swaption,
};
#[allow(unused_imports)]
use ql_rates::{
    HuntKennedyMarkovFunctionalReport, HuntKennedyMarkovFunctionalSpec,
    evaluate_hunt_kennedy_markov_functional_range_accrual,
};
use ql_rates::{
    TurfusWillemsRfrInArrearsSmileSpec, evaluate_turfus_willems_rfr_in_arrears_caplet_smile,
};
use ql_time::{BitsetCalendar, BusinessDayConvention, Date, DayCounter, Period};
#[allow(unused_imports)]
use ql_vol::{
    BergomiForwardStartNapoleonVolReport, BergomiForwardStartNapoleonVolSpec,
    evaluate_bergomi_forward_start_smile_and_napoleon_convexity,
};
use ql_vol::{
    BourgeyGatheralBergomiGuyonSmileSpec, evaluate_bourgey_gatheral_bergomi_guyon_expansion,
};
use ql_vol::{
    CalibratedSabrSmile, FreeBoundarySabrParams, SmileStrikeQuote, ThirdOrderSpeedVolExpansionSpec,
    evaluate_free_boundary_sabr, evaluate_third_order_speed_vol_expansion, implied_black_vol,
};
#[allow(unused_imports)]
use ql_vol::{
    GatheralDemeterfiCorridorVarianceSkewReport, GatheralDemeterfiCorridorVarianceSkewSpec,
    evaluate_gatheral_demeterfi_corridor_variance_and_skew,
};
#[allow(unused_imports)]
use ql_vol::{
    MingoneGlobalArbitrageFreeEssviReport, MingoneGlobalArbitrageFreeEssviSpec,
    evaluate_mingone_global_arbitrage_free_essvi,
};
#[allow(unused_imports)]
use ql_xva::{
    AlbaneseCrepeyMporBalanceSheetReport, AlbaneseCrepeyMporBalanceSheetSpec,
    compute_albanese_crepey_mpor_and_balance_sheet_xva,
};
use ql_xva::{
    AlbaneseCrepeyZeronBalanceSheetXvaSpec, evaluate_albanese_crepey_zeron_balance_sheet_xva,
};
use ql_xva::{
    BaCvaCounterpartySpec, BaCvaPortfolioSpec, BaCvaSingleNameHedge, DynamicSimmMvaSpec,
    MultiCurrencyCsaSpec, NettingSetSwap, SaCcrAssetClass, SaCcrNettingSet, SaCcrTrade, SimmBucket,
    TurlakovWwrSpec, compute_ba_cva_capital, compute_dynamic_simm_mva_and_spike, compute_frtb_sa,
    compute_isda_simm, compute_isda_simm_6_risk_classes, compute_multicurrency_ctd_colva,
    compute_sa_ccr_ead, compute_turlakov_wwr_cva_fva, compute_xva, expected_positive_exposure,
    potential_future_exposure, simulate_hw1f_exposure_cube,
};
#[allow(unused_imports)]
use ql_xva::{
    BcbsD424OutputFloorAndSmaOpRiskReport, BcbsD424OutputFloorAndSmaOpRiskSpec,
    evaluate_bcbs_d424_output_floor_and_sma_op_risk,
};
#[allow(unused_imports)]
use ql_xva::{
    GreenFrtbCvaHedgeCapitalReport, GreenFrtbCvaHedgeCapitalSpec,
    evaluate_green_frtb_cva_hedge_relief_and_kva,
};
use std::fmt::Write as FmtWrite;
use std::sync::Arc;
use std::time::Instant;

const BRAILLE_DOTS: [[u32; 2]; 4] = [[0x01, 0x08], [0x02, 0x10], [0x04, 0x20], [0x40, 0x80]];

/// Logical workstation tabs.  The dashboard remains headless and printable,
/// while callers can use this stable enumeration to select a domain view.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum TabId {
    Market,
    Rates,
    Credit,
    Inflation,
    Fx,
    Equities,
    Volatility,
    Xva,
    Portfolio,
    Exotics,
    Commodities,
    Levy,
    FixedIncome,
    MonteCarlo,
    Analytics,
    Blotter,
    AmericanHybrid,
    BlRiskParity,
    BermudanVolIft,
    SaCcrSwingQp,
    SlvRoughWorstOf,
    FrtbCvarExecution,
    MonotoneConvexVixParisianRobust,
    DoubleHestonTarfCocoHrp,
    JyInflationTbaBergomiNco,
    CheyetteCmoQrhBaseCorrCdar,
    AswMf1fAndreasenHugeBessEp,
    MtmXccsColvaZabrPdvWeatherEnb,
    BurnoutCmcdsQuinticAfvCsSchur,
    SinHouMbRBergomiChoiAfsExpectile,
    MonotoneConvexClarkSaccrRvolCartea,
    KrdLpiCocoBergomi2fCgmyGueantMeucci,
    RkhsFmmVvLsvJltCtdQrhParisianAidNco,
    MbsIoPoRmvMoXvaDilutionSigGarchPot,
    KalmanRtsVecmHjSdfIfrs9CptHasbrouck,
    MfRatesWildcardCtdPivotTarfBergomi,
    CochraneFtplG2TwistEssviOffshoreCb,
    BcbsIrrbbTurfusRfrBaranMtmBourgeyBg,
    VeronesiConvLandoGenGatheralCorridor,
    RomerTaylorDamodaranCrpChristoffersenVar,
    FabozziKrdCmoPacAsquithCbAsc470Colva,
    CheyetteCapletAdbNdirsZabrQrhCmcdsAdi,
    RebonatoCmsLelandToftAscotPoolRollRie,
    TuckmanRepoBkFvaDecsIoPoKylePinGjrLppl,
    HaganWestCtdSwitchAfvCbKthDefaultSaCva,
    BrigoCsaCtdCallableOasJltCorridorLsm,
    BrigoNssKrdTipsFloorRmvTriangleStickyLv,
    HaganWestTensionPacCmoFamaMacBethAsviBgk,
    BrigoXccyMtmSpecPoolMarshallOlkinVixVvix,
    BrigoOisFallbackCollaredFrnCdsFepSparkPrdc,
    AndersenCmsSpecialRepoCdsUpfrontSsrKvaCrack,
    DarbyshireStepOisCbSoftCallCdxSkewMalzFly,
    HenrardTenorBasisMsrRecaptureCdoBaseCorrKiko,
    BrigoDiffSwapCallableMuniBkmSkewBreakCva,
    AndersenSofrFutHayreMbsPedersenCdxDupireLv,
    LittermanScheinkmanPcaFlyTfSplitCbDoubleTCdoLeeWings,
    PiterbargCollateralOptTuckmanRolldownMertonKmvBartlettSabr,
    GswSvenssonCurveFabozziMbsEsfVasicekGordyIrbCarrWuVol,
    KimWrightShadowBrennanSchwartzCbDuffieSingletonRmvFukasawaGamma,
    JpsUnspannedMacroKalotayRefundingHullWhiteNtdDermanSticky,
    LyashenkoFmmCapletTuckmanCpiSeasonalityBrigoAlfonsiCirFouqueSv,
    JszCanonicalGdtsmTsiveriotisFernandesCbLaurentGregoryCdoBroadieJainVarswap,
    BianchettiTwoCurveQuantoBdtCallableBondSchonbucherCdsSwaptionBatesSvj,
    AndersenPiterbargCmsSpreadFabozziZBondCmoDuffieGarleanuAffineCdoCarrLee,
    PiterbargCtdCollateralTuckmanMakeWholeBondBrigoMoriniCmcdsHaganShiftedSabr,
    MercurioFvaMultiCurveSwapIngersollContingentCallCbHullWhiteImpliedCopulaCdoGatheralJacquierSviJw,
    AmetranoBianchettiTenorBasisFabozziHayreIoPoStripHullBrigoWwrCdsBossensVannaVolgaBarrier,
    SmithWilsonUfrKalotayWilliamsFabozziCallableSchonbucherClaytonCdoBuehlerVarianceCurve,
    McCullochTaxCurveAndersenBuffumCbOkaneRecoveryLockCdsDermanKaniImpliedTree,
    HaganWestMonotoneConvexFabozziEffectiveDurationOasFryeAltmanDownturnLgdGatheralJacquierSviArbitrage,
    MercurioRussoIborOisBasisJarrowYildirimYoyInflationDuffieLandoHybridCdsBayerFrizGatheralRoughSkew,
    FilipovicLarssonTrolleLinearRationalDuffieKrishnamurthySpecialRepoDuffieSaitaWangStochasticDefaultCarrMadanNeubergerLogContract,
    LyashenkoMercurioRfrInArrearsPennacchiSundaresanCocoAt1PedersenBrigoCdsIndexFepBroadieJainSeppSvjVolswap,
    CieslakPovalaTrendInflationCycleGabaixKrishnamurthyVigneronMbsOasErraisGieseckeGoldbergHawkesCreditMartinBondarenkoSvixLowerBound,
    BianchettiMercurioMultiCurveFraQuantoAcharyaAmihudBharathLiquidityStressBrigoAlfonsiCorrelatedCirCdsFukasawaAsymptoticSkewVarswap,
    HansonSteinGreenwoodVayanosDurationSupplyFleckensteinLongstaffLustigTipsPuzzleCollinDufresneGoldsteinYangSuperSeniorCdoTodorovOptionSpotVolatility,
    KlinglerSyrstadSofrRepoSpikeLongstaffMithalNeisCdsBondBasisPanSingletonSovereignCdsMedvedevScailletShortMaturitySmile,
}
impl TabId {
    pub const TAB_COUNT: usize = 77;
    pub const ALL: [Self; 77] = [
        Self::Market,
        Self::Rates,
        Self::Credit,
        Self::Inflation,
        Self::Fx,
        Self::Equities,
        Self::Volatility,
        Self::Xva,
        Self::Portfolio,
        Self::Exotics,
        Self::Commodities,
        Self::Levy,
        Self::FixedIncome,
        Self::MonteCarlo,
        Self::Analytics,
        Self::Blotter,
        Self::AmericanHybrid,
        Self::BlRiskParity,
        Self::BermudanVolIft,
        Self::SaCcrSwingQp,
        Self::SlvRoughWorstOf,
        Self::FrtbCvarExecution,
        Self::MonotoneConvexVixParisianRobust,
        Self::DoubleHestonTarfCocoHrp,
        Self::JyInflationTbaBergomiNco,
        Self::CheyetteCmoQrhBaseCorrCdar,
        Self::AswMf1fAndreasenHugeBessEp,
        Self::MtmXccsColvaZabrPdvWeatherEnb,
        Self::BurnoutCmcdsQuinticAfvCsSchur,
        Self::SinHouMbRBergomiChoiAfsExpectile,
        Self::MonotoneConvexClarkSaccrRvolCartea,
        Self::KrdLpiCocoBergomi2fCgmyGueantMeucci,
        Self::RkhsFmmVvLsvJltCtdQrhParisianAidNco,
        Self::MbsIoPoRmvMoXvaDilutionSigGarchPot,
        Self::KalmanRtsVecmHjSdfIfrs9CptHasbrouck,
        Self::MfRatesWildcardCtdPivotTarfBergomi,
        Self::CochraneFtplG2TwistEssviOffshoreCb,
        Self::BcbsIrrbbTurfusRfrBaranMtmBourgeyBg,
        Self::VeronesiConvLandoGenGatheralCorridor,
        Self::RomerTaylorDamodaranCrpChristoffersenVar,
        Self::FabozziKrdCmoPacAsquithCbAsc470Colva,
        Self::CheyetteCapletAdbNdirsZabrQrhCmcdsAdi,
        Self::RebonatoCmsLelandToftAscotPoolRollRie,
        Self::TuckmanRepoBkFvaDecsIoPoKylePinGjrLppl,
        Self::HaganWestCtdSwitchAfvCbKthDefaultSaCva,
        Self::BrigoCsaCtdCallableOasJltCorridorLsm,
        Self::BrigoNssKrdTipsFloorRmvTriangleStickyLv,
        Self::HaganWestTensionPacCmoFamaMacBethAsviBgk,
        Self::BrigoXccyMtmSpecPoolMarshallOlkinVixVvix,
        Self::BrigoOisFallbackCollaredFrnCdsFepSparkPrdc,
        Self::AndersenCmsSpecialRepoCdsUpfrontSsrKvaCrack,
        Self::DarbyshireStepOisCbSoftCallCdxSkewMalzFly,
        Self::HenrardTenorBasisMsrRecaptureCdoBaseCorrKiko,
        Self::BrigoDiffSwapCallableMuniBkmSkewBreakCva,
        Self::AndersenSofrFutHayreMbsPedersenCdxDupireLv,
        Self::LittermanScheinkmanPcaFlyTfSplitCbDoubleTCdoLeeWings,
        Self::PiterbargCollateralOptTuckmanRolldownMertonKmvBartlettSabr,
        Self::GswSvenssonCurveFabozziMbsEsfVasicekGordyIrbCarrWuVol,
        Self::KimWrightShadowBrennanSchwartzCbDuffieSingletonRmvFukasawaGamma,
        Self::JpsUnspannedMacroKalotayRefundingHullWhiteNtdDermanSticky,
        Self::LyashenkoFmmCapletTuckmanCpiSeasonalityBrigoAlfonsiCirFouqueSv,
        Self::JszCanonicalGdtsmTsiveriotisFernandesCbLaurentGregoryCdoBroadieJainVarswap,
        Self::BianchettiTwoCurveQuantoBdtCallableBondSchonbucherCdsSwaptionBatesSvj,
        Self::AndersenPiterbargCmsSpreadFabozziZBondCmoDuffieGarleanuAffineCdoCarrLee,
        Self::PiterbargCtdCollateralTuckmanMakeWholeBondBrigoMoriniCmcdsHaganShiftedSabr,
        Self::MercurioFvaMultiCurveSwapIngersollContingentCallCbHullWhiteImpliedCopulaCdoGatheralJacquierSviJw,
        Self::AmetranoBianchettiTenorBasisFabozziHayreIoPoStripHullBrigoWwrCdsBossensVannaVolgaBarrier,
        Self::SmithWilsonUfrKalotayWilliamsFabozziCallableSchonbucherClaytonCdoBuehlerVarianceCurve,
        Self::McCullochTaxCurveAndersenBuffumCbOkaneRecoveryLockCdsDermanKaniImpliedTree,
        Self::HaganWestMonotoneConvexFabozziEffectiveDurationOasFryeAltmanDownturnLgdGatheralJacquierSviArbitrage,
        Self::MercurioRussoIborOisBasisJarrowYildirimYoyInflationDuffieLandoHybridCdsBayerFrizGatheralRoughSkew,
        Self::FilipovicLarssonTrolleLinearRationalDuffieKrishnamurthySpecialRepoDuffieSaitaWangStochasticDefaultCarrMadanNeubergerLogContract,
        Self::LyashenkoMercurioRfrInArrearsPennacchiSundaresanCocoAt1PedersenBrigoCdsIndexFepBroadieJainSeppSvjVolswap,
        Self::CieslakPovalaTrendInflationCycleGabaixKrishnamurthyVigneronMbsOasErraisGieseckeGoldbergHawkesCreditMartinBondarenkoSvixLowerBound,
        Self::BianchettiMercurioMultiCurveFraQuantoAcharyaAmihudBharathLiquidityStressBrigoAlfonsiCorrelatedCirCdsFukasawaAsymptoticSkewVarswap,
        Self::HansonSteinGreenwoodVayanosDurationSupplyFleckensteinLongstaffLustigTipsPuzzleCollinDufresneGoldsteinYangSuperSeniorCdoTodorovOptionSpotVolatility,
        Self::KlinglerSyrstadSofrRepoSpikeLongstaffMithalNeisCdsBondBasisPanSingletonSovereignCdsMedvedevScailletShortMaturitySmile,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Self::Market => "01 Market",
            Self::Rates => "02 Rates",
            Self::Credit => "03 Credit",
            Self::Inflation => "04 Inflation",
            Self::Fx => "05 FX",
            Self::Equities => "06 Equities",
            Self::Volatility => "07 Volatility",
            Self::Xva => "08 XVA",
            Self::Portfolio => "09 Portfolio",
            Self::Exotics => "10 Exotics",
            Self::Commodities => "11 Commodities",
            Self::Levy => "12 Levy",
            Self::FixedIncome => "13 Fixed Income",
            Self::MonteCarlo => "14 Monte Carlo",
            Self::Analytics => "15 Analytics",
            Self::Blotter => "16 Blotter",
            Self::AmericanHybrid => "17 American / Hybrid",
            Self::BlRiskParity => "18 BL / Risk Parity",
            Self::BermudanVolIft => "19 Bermudan / Vol IFT",
            Self::SaCcrSwingQp => "20 SA-CCR / Swing / QP",
            Self::SlvRoughWorstOf => "21 SLV / rBergomi / WoF",
            Self::FrtbCvarExecution => "22 FRTB / CVaR / Exec",
            Self::MonotoneConvexVixParisianRobust => "23 Monotone Convex / VIX / Parisian / Robust",
            Self::DoubleHestonTarfCocoHrp => "24 Double Heston / TARF / CoCo / HRP",
            Self::JyInflationTbaBergomiNco => "25 JY-Inflation / TBA-Roll / Bergomi2F / NCO",
            Self::CheyetteCmoQrhBaseCorrCdar => "26 Cheyette / CMO / QRH / BaseCorr / EVaR",
            Self::AswMf1fAndreasenHugeBessEp => "27 ASW / MF1F / AH-LV / BESS / Idzorek-DP",
            Self::MtmXccsColvaZabrPdvWeatherEnb => {
                "28 MtM-XCCS / CTD-COLVA / ZABR-PDV / Weather / ENB"
            }
            Self::BurnoutCmcdsQuinticAfvCsSchur => {
                "29 Burnout-MBS / CMCDS / Quintic-OU / AFV-CB / CS-Swaption / Schur-HRP"
            }
            Self::SinHouMbRBergomiChoiAfsExpectile => {
                "30 Sinh-OU / BK-FRN / Morini-Brigo / WWM-XVA / mrBergomi / Choi / AFS-LOB / EXVaR"
            }
            Self::MonotoneConvexClarkSaccrRvolCartea => {
                "31 HW-MC / Clark-FX / SA-CCR-KVA / BFG-rVol / Cartea-Alpha / Paleologo"
            }
            Self::KrdLpiCocoBergomi2fCgmyGueantMeucci => {
                "32 KRD / LPI / CoCo / Bergomi / CGMY / AAD"
            }
            Self::RkhsFmmVvLsvJltCtdQrhParisianAidNco => "33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML",
            Self::MbsIoPoRmvMoXvaDilutionSigGarchPot => {
                "34 MBS-IO-PO / RMV / MO-XVA / Dilution / Sig / GARCH / POT / RIE / PIN / GP / Brinson / IERP"
            }
            Self::KalmanRtsVecmHjSdfIfrs9CptHasbrouck => {
                "35 Kalman-RTS / VECM / HJ-SDF / IFRS9 / CPT / Hasbrouck / FV / Lando / LPPL / CB-Arb / Causal / Malliavin"
            }
            Self::MfRatesWildcardCtdPivotTarfBergomi => {
                "36 MF-Rates / Wildcard / TARF / Bergomi / Schönbucher / FRTB-CVA / Lookback / VTA / FLAM / ACD / HLZ / MLMC"
            }
            Self::CochraneFtplG2TwistEssviOffshoreCb => {
                "37 FTPL / G2++ / eSSVI / Offshore-CB / CoCo / MPOR / Double-Heston / Curran / Warin / Roll-Corwin / Busseti-Kelly / Diff-QP"
            }
            Self::BcbsIrrbbTurfusRfrBaranMtmBourgeyBg => {
                "38 IRRBB / Turfus-RFR / Baran-MtM / Bourgey-BG / FEDS-PONV / Zeron-SIMM / Jacquier-rBergomi / Choi-Spread / BESS-SDP / Kratz-Dark / Föllmer / Diff-PCA"
            }
            Self::VeronesiConvLandoGenGatheralCorridor => {
                "39 Veronesi / Lando-Log / Corridor-Var / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF-AP / IC-Decay / Glasserman-IS"
            }
            Self::RomerTaylorDamodaranCrpChristoffersenVar => {
                "40 Taylor-Curve / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect-CGO / Regret-Match / Beneish-HKAS32 / Campisi-FI / Heckman-GMM / Elastic-Net / Reset-CB / FSB-BailIn"
            }
            Self::FabozziKrdCmoPacAsquithCbAsc470Colva => {
                "41 KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI"
            }
            Self::CheyetteCapletAdbNdirsZabrQrhCmcdsAdi => {
                "42 Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan"
            }
            Self::RebonatoCmsLelandToftAscotPoolRollRie => {
                "43 Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu"
            }
            Self::TuckmanRepoBkFvaDecsIoPoKylePinGjrLppl => {
                "44 Tuckman-Repo / BK-FVA / DECS / MBS-IO-PO / Kyle-PIN / GP-Aim / GJR-DCC / LPPL / Spectral-Risk / Brinson-Carino / Prelec-CPT / Path-Sig"
            }
            Self::HaganWestCtdSwitchAfvCbKthDefaultSaCva => {
                "45 Hagan-West-MC / CTD-Switch / AFV-CB / Kth-Default / SA-CVA-Euler / GS-Kirk-Spread / FX-Vanna-Volga / Napoleon-Cliquet / Transfer-Coef / Meucci-Entropy / KRR-LOOCV / McKinsey-EVA"
            }
            Self::BrigoCsaCtdCallableOasJltCorridorLsm => {
                "46 Multi-CSA-CTD / Callable-OAS / JLT-Migration / Corridor-VarSwap / LSM-Rogers-Dual / Cartea-LOB-MM / Lo-MacKinlay-VR / Lopez-DSR / Roncalli-ERC / Ilmanen-VMC / Penman-ReOI / McCall-Bellman"
            }
            Self::BrigoNssKrdTipsFloorRmvTriangleStickyLv => {
                "47 NSS-KRD / TIPS-Floor-BEI / RMV-RT-RFV / Sticky-LocalVol / SIMM-CCP-MVA / Schwartz-Smith-KF / FX-DNT-VV / Timer-VarClock / BL-Woodbury / AC-Sinh-Liq / ElasticNet-Group / NK-DSGE-REE"
            }
            Self::HaganWestTensionPacCmoFamaMacBethAsviBgk => {
                "48 Tension-Curve / Duffie-Lando / PAC-CMO / Fama-MacBeth / Gatheral-aSVI / BGK-Barrier / Cross-Impact / Karnosky-Singer / Merton-JD / Campbell-Viceira / Ohlson-LID / Revolve-BSB"
            }
            Self::BrigoXccyMtmSpecPoolMarshallOlkinVixVvix => {
                "49 XCCY-MtM / Spec-Pool-IO-PO / Marshall-Olkin / VIX-VVIX / WWR-CVA / Gibson-Schwartz / FX-FVA / Chooser / FLAM-TC / ES-Euler / James-Stein-SURE / Campbell-Cochrane"
            }
            Self::BrigoOisFallbackCollaredFrnCdsFepSparkPrdc => {
                "50 OIS-Fallback / Collared-FRN / CDS-Swaption-FEP / Bergomi-2F / SA-CCR-EAD / Clean-Spark / FX-PRDC / Geske-Compound / Alpha-Decay-IR / Meucci-ENB / Johansen-VECM / Bansal-Yaron"
            }
            Self::AndersenCmsSpecialRepoCdsUpfrontSsrKvaCrack => {
                "51 CMS-Replication / Special-Repo-CTD / ISDA-CDS-Upfront / SSR-Tanaka / SA-CCR-KVA / 3-Asset-Crack / PRDC-Cap / ESO-Dilution / Spinu-ERC / Rockafellar-CVaR / PLS-Krylov / FOMC-Hike"
            }
            Self::DarbyshireStepOisCbSoftCallCdxSkewMalzFly => {
                "52 Step-OIS-Turn / CB-SoftCall-MW / CDX-Skew-Basis / VarSwap-Trunc / FVA-DVA-Overlap / Working-Storage / FX-Malz-Fly / GSG-Lookback / Ledoit-Wolf-CC / Almgren-CARA / ESTAR-Band / ACM-CP-TP"
            }
            Self::HenrardTenorBasisMsrRecaptureCdoBaseCorrKiko => {
                "53 Tenor-Basis-Mult / MSR-Recapture / CDO-Base-Corr / Fwd-Start-Smile / CSA-MTA-CVA / Gas-Swing / FX-Window-KIKO / Stulz-Rainbow / Jagannathan-Ma / Obizhaeva-Wang / Newey-West-HAC / Gabaix-Disaster"
            }
            Self::BrigoDiffSwapCallableMuniBkmSkewBreakCva => {
                "54 Diff-Swap-Quanto / Callable-Muni-YTW / Step-Up-Bond / BKM-SKEW / Break-Clause-CVA / Power-Skew / FX-Fader / Longstaff-Extendible / Kan-Zhou-3Fund / Sharpe-RBSA / Diebold-Mariano-HLN / Laubach-Williams-r*"
            }
            Self::AndersenSofrFutHayreMbsPedersenCdxDupireLv => {
                "55 SOFR-Fut-FRA / Hayre-MBS-SCurve / CDX-Option-FEP / Dupire-TotalVar-LV / SIMM-MVA-Funding / RSS-Stockout / CIP-Basis-Wedge / Carr-Chou-Symmetry / Pastor-Stambaugh / Treynor-Black / Engle-Granger-ECM / Epstein-Zin-SDF"
            }
            Self::LittermanScheinkmanPcaFlyTfSplitCbDoubleTCdoLeeWings => {
                "56 PCA-Fly-Carry / TF-Split-CB / Double-t-CDO / Roger-Lee-Wings / CDS-Contagion-CVA / Schwartz-OU-Fut / FX-DCD-Yield / Passport-Opt / GRS-Factor-FTest / Sqrt-Metaorder / BNS-HAR-CJ / Lettau-Ludvigson-cay"
            }
            Self::PiterbargCollateralOptTuckmanRolldownMertonKmvBartlettSabr => {
                "57 Piterbarg-Coll-Opt / Tuckman-Rolldown / Merton-KMV-DD / Bartlett-SABR / Pykhtin-MPoR-Spike / Deng-FTR-Opt / Wystup-FX-Digital / Conze-Fixed-Lookback / Tasche-Hot-Spots-ES / Frazzini-BAB / Pesaran-Timmermann / Lucas-Mehra-Prescott"
            }
            Self::GswSvenssonCurveFabozziMbsEsfVasicekGordyIrbCarrWuVol => {
                "58 GSW-Svensson-Curve / Fabozzi-MBS-ESF / Vasicek-Gordy-IRB-GA / Carr-Wu-Vol-Drift / Gregory-Marginal-CVA / Working-Brennan-Storage / Wystup-Forward-Extra / Rubinstein-Exchange-Opt / HLZ-Factor-FDR / Holden-Subrahmanyam-Kyle / Bollerslev-GARCH-Kurtosis / Campbell-Shiller-CAPE"
            }
            Self::KimWrightShadowBrennanSchwartzCbDuffieSingletonRmvFukasawaGamma => {
                "59 Kim-Wright-Shadow / Brennan-Schwartz-CB / Duffie-Singleton-RMV / Fukasawa-Gamma-Swap / Burgard-Kjaer-FVA / Carmona-Durrleman-Spark / Wystup-FX-TARF / BGK-Discrete-Lookback / Novy-Marx-FF5 / Acerbi-Spectral-ES / Corsi-HAR-RV / BLL-News-Noise"
            }
            Self::JpsUnspannedMacroKalotayRefundingHullWhiteNtdDermanSticky => {
                "60 JPS-Unspanned-Macro / Kalotay-Refunding / Hull-White-NtD / Derman-Sticky-Regimes / Hull-White-WWR-CVA / Gibson-Schwartz-Samuelson / Reiner-Quanto-FX / Turnbull-Wakeman-Asian / AMP-Value-Momentum / Garleanu-Pedersen-Aim / Bai-Ng-ICp2 / BGG-Accelerator"
            }
            Self::LyashenkoFmmCapletTuckmanCpiSeasonalityBrigoAlfonsiCirFouqueSv => {
                "61 Lyashenko-FMM-Caplet / Tuckman-CPI-Seasonality / Brigo-Alfonsi-CIR-CDS / Fouque-Multiscale-SV / Gregory-Bilateral-FtD / Eydeland-VPP-Heat-Rate / Wystup-FX-Corridor / Rubinstein-Reiner-Barrier / LNS-GLS-R2 / Bellini-Expectile / BNS-Quarticity / Kiyotaki-Moore-Cycle"
            }
            Self::JszCanonicalGdtsmTsiveriotisFernandesCbLaurentGregoryCdoBroadieJainVarswap => {
                "62 JSZ-Canonical-GDTSM / Tsiveriotis-Fernandes-CB / Laurent-Gregory-CDO / Broadie-Jain-VarSwap / Albanese-SA-CCR-KVA / RSS-Equilibrium-Storage / Du-Tepper-Verdelhan-CIP / Wilmott-Cliquet / BSV-Parametric-Portfolio / Ahmadi-Javid-EVaR / Engle-Kroner-BEKK / Gabaix-Rare-Disaster"
            }
            Self::BianchettiTwoCurveQuantoBdtCallableBondSchonbucherCdsSwaptionBatesSvj => {
                "63 Bianchetti-Two-Curve-Quanto / BDT-Callable-Bond / Schonbucher-CDS-Swaption / Bates-SVJ-Smile / Brigo-Pallavicini-CCP-XVA / DMT-Cointegrated-Crack / Piterbarg-PRDC-TARN / BGK-Discrete-Barrier / Connor-Korajczyk-APT / AFS-Transient-Impact / Hamilton-Markov-Switching / Bansal-Yaron-LRR"
            }
            Self::AndersenPiterbargCmsSpreadFabozziZBondCmoDuffieGarleanuAffineCdoCarrLee => {
                "64 Andersen-Piterbarg-CMS-Spread / Fabozzi-Z-Bond-CMO / Duffie-Garleanu-Affine-CDO / Carr-Lee-Vol-Swap / Green-Leverage-Ratio-KVA / Cortazar-Schwartz-3F-Oil / Fujii-Collateralized-XCCS / Pechtl-Haug-Time-Switch / Boudt-Cornish-Fisher-MVaR / Rockafellar-Quadrangle / Hasbrouck-Info-Share / He-Krishnamurthy-SDF"
            }
            Self::PiterbargCtdCollateralTuckmanMakeWholeBondBrigoMoriniCmcdsHaganShiftedSabr => {
                "65 Piterbarg-CTD-Collateral / Tuckman-Make-Whole-Bond / Brigo-Morini-CMCDS / Hagan-Shifted-SABR / Green-Kenyon-TVA / Eydeland-Interconnector-Spread / Wystup-Clark-FX-Strangle / Kunitomo-Ikeda-Curved-Barrier / CLMX-Variance-Decomp / Kyle-Multi-Period-Lambda / Hansen-Lunde-Overnight-RV / Goldstein-Pauzner-Bank-Run"
            }
            Self::MercurioFvaMultiCurveSwapIngersollContingentCallCbHullWhiteImpliedCopulaCdoGatheralJacquierSviJw => {
                "66 Mercurio-FVA-Swap / Ingersoll-Contingent-Call-CB / Hull-White-Implied-Copula-CDO / Gatheral-Jacquier-SVI-JW / Burgard-Kjaer-Semi-Rep-XVA / Benth-Seasonal-HDD-CDD-Swap / Baba-Packer-FX-Swap-Basis / Conze-Viswanathan-Lookback / He-Litterman-BL-Posterior / GHS-Importance-Sampling-CVaR / Engle-Granger-Johansen-VECM / Gabaix-Rare-Disasters"
            }
            Self::AmetranoBianchettiTenorBasisFabozziHayreIoPoStripHullBrigoWwrCdsBossensVannaVolgaBarrier => {
                "67 Ametrano-Bianchetti-Tenor-Basis / Fabozzi-Hayre-IO-PO-Strip / Hull-Brigo-WWR-CDS / Bossens-Vanna-Volga-Barrier / Gregory-SIMM-MVA / Carmona-Clean-Dark-Spread / Reiner-Haug-Quanto-Call / Hart-Ross-Soft-Barrier / Pastor-Stambaugh-Alpha-Shrinkage / Obizhaeva-Wang-Impulse-LOB / Diebold-Mariano-West-HLN / Piazzesi-Schneider-Nominal-TP"
            }
            Self::SmithWilsonUfrKalotayWilliamsFabozziCallableSchonbucherClaytonCdoBuehlerVarianceCurve => {
                "68 Smith-Wilson-UFR / Kalotay-Williams-Fabozzi-Callable / Schonbucher-Clayton-CDO / Buehler-Variance-Curve / Pykhtin-Rosen-Marginal-CVA / Hilliard-Reis-Jump-Futures / Gabaix-Maggiori-FX-Carry / Cheuk-Vorst-Shout-Option / Grossman-Zhou-Drawdown-Kelly / GLS-Tasche-Cornish-Fisher-RAROC / Andrews-Ploberger-Wald-Break / Jermann-Quadrini-Financial-Shock"
            }
            Self::McCullochTaxCurveAndersenBuffumCbOkaneRecoveryLockCdsDermanKaniImpliedTree => {
                "69 McCulloch-Tax-Curve / Andersen-Buffum-CB / OKane-Recovery-Lock-CDS / Derman-Kani-Implied-Tree / Kenyon-Green-CSA-LVA / Boogert-de-Jong-Gas-Ratchet / Garman-Kohlhagen-BAW-FX / Hui-Haug-Binary-Barrier / Kan-Robotti-HJ-Distance / Cartea-Jaimungal-Order-Flow / Ait-Sahalia-Jacod-Jump-Test / Wachter-Disaster-Volatility"
            }
            Self::HaganWestMonotoneConvexFabozziEffectiveDurationOasFryeAltmanDownturnLgdGatheralJacquierSviArbitrage => {
                "70 Hagan-West-Monotone-Convex / Fabozzi-Effective-Duration-OAS / Frye-Altman-Downturn-LGD / Gatheral-Jacquier-SVI-Arbitrage / Gregory-Brigo-FtD-Bilateral-CVA / Brennan-Schwartz-Mine-Hysteresis / Wystup-FX-Fader / Broadie-Detemple-Capped-Call / Meucci-Entropy-Pooling / Cont-Deguest-Scandolo-Robust-Risk / Pesaran-Shin-Smith-ARDL / Christiano-Eichenbaum-Evans-NK"
            }
            Self::MercurioRussoIborOisBasisJarrowYildirimYoyInflationDuffieLandoHybridCdsBayerFrizGatheralRoughSkew => {
                "71 Mercurio-Russo-IBOR-OIS-Basis / Jarrow-Yildirim-YoY-Inflation / Duffie-Lando-Hybrid-CDS / Bayer-Friz-Gatheral-Rough-Skew / Andersen-Pykhtin-Sokol-MPoR / Bjerksund-Stensland-321-Crack / Carr-Chou-Symmetry-Barrier / Chesney-Jeanblanc-Yor-Parisian / DeMiguel-Brodie-Norm-Shrinkage / Gueant-Lehalle-Tapia-MM / Hansen-White-SPA-Test / Brunnermeier-Sannikov-Paradox"
            }
            Self::FilipovicLarssonTrolleLinearRationalDuffieKrishnamurthySpecialRepoDuffieSaitaWangStochasticDefaultCarrMadanNeubergerLogContract => {
                "72 Filipovic-Larsson-Trolle-LRTSM / Duffie-Krishnamurthy-Repo / Duffie-Saita-Wang-Default / Carr-Madan-Neuberger-LogContract / Hull-White-Burgard-Kjaer-FVA / Casassus-Collin-Dufresne-CY / Lustig-Roussanov-Verdelhan-FX / Carr-Lee-Corridor-VarSwap / Garlappi-Uppal-Wang-MultiPrior / McNeil-Frey-EVT-GARCH / Patton-SJC-Copula / Pastor-Veronesi-Political-Risk"
            }
            Self::LyashenkoMercurioRfrInArrearsPennacchiSundaresanCocoAt1PedersenBrigoCdsIndexFepBroadieJainSeppSvjVolswap => {
                "73 Lyashenko-Mercurio-RFR / Pennacchi-Sundaresan-CoCo / Pedersen-Brigo-CDX-FEP / Broadie-Jain-Sepp-SVJ-VolSwap / Brigo-Capponi-Pallavicini-WWR / Routledge-Seppi-Spatt-Storage / Piterbarg-PRDC-Hybrid / Broadie-Glasserman-Kou-Barrier / Brandt-Santa-Clara-Valkanov-PPP / Gatheral-Schied-Slynko-Impact / Corsi-BPQ-HARQ / Lettau-Ludvigson-CAY"
            }
            Self::CieslakPovalaTrendInflationCycleGabaixKrishnamurthyVigneronMbsOasErraisGieseckeGoldbergHawkesCreditMartinBondarenkoSvixLowerBound => {
                "74 Cieslak-Povala-Cycle / Gabaix-Krishnamurthy-Vigneron-MBS / Errais-Giesecke-Goldberg-Hawkes / Martin-Bondarenko-SVIX / Crepey-Bichuch-Capponi-BSDE / Trolle-Schwartz-USV / Ivashina-Scharfstein-Stein-CIP / Wilmott-Windcliff-Cliquet / Ledoit-Wolf-Nonlinear-Shrinkage / Adrian-Brunnermeier-CoVaR / Bollerslev-Tauchen-Zhou-VRP / Campbell-Vuolteenaho-Beta"
            }
            Self::BianchettiMercurioMultiCurveFraQuantoAcharyaAmihudBharathLiquidityStressBrigoAlfonsiCorrelatedCirCdsFukasawaAsymptoticSkewVarswap => {
                "75 Bianchetti-Mercurio-FRA-Quanto / Acharya-Amihud-Bharath-Liq / Brigo-Alfonsi-SSRD-CDS / Fukasawa-Skew-VarSwap / Andersen-Duffie-Song-FVA / Pirrong-Inventory-Spread / Lipton-Castagna-FX-DNT / Linetsky-Step-Barrier / Koijen-Yogo-Demand-System / Kyle-Obizhaeva-Invariance / Gu-Kelly-Xiu-OOS-R2 / Drechsler-Savov-Schnabl-Deposit"
            }
            Self::HansonSteinGreenwoodVayanosDurationSupplyFleckensteinLongstaffLustigTipsPuzzleCollinDufresneGoldsteinYangSuperSeniorCdoTodorovOptionSpotVolatility => {
                "76 Hanson-Stein-Duration / Fleckenstein-Longstaff-Lustig-TIPS / Collin-Dufresne-Goldstein-Yang-CDO / Todorov-Spot-Vol / Albanese-Caenazzo-Crepey-KVA / Acharya-Lochstoer-Ramadorai-Hedging / Brunnermeier-Nagel-Pedersen-FX / Broadie-Yamamoto-FGT-Asian / Kelly-Pruitt-Su-IPCA / Brownlees-Engle-SRISK / ABDL-ARFIMA-RV / Farhi-Gabaix-Disaster-FX"
            }
            Self::KlinglerSyrstadSofrRepoSpikeLongstaffMithalNeisCdsBondBasisPanSingletonSovereignCdsMedvedevScailletShortMaturitySmile => {
                "77 Klingler-Syrstad-SOFR-Spike / Longstaff-Mithal-Neis-CDS-Bond / Pan-Singleton-Sovereign-CDS / Medvedev-Scaillet-Short-Smile / Duffie-Zhu-CCP-Netting / Gorton-Hayashi-Rouwenhorst-Inventory / Della-Corte-Ramadorai-Sarno-FX-VRP / Pelsser-Laplace-Double-Barrier / Moreira-Muir-Vol-Managed / Budish-Cramton-Shim-FBA / AMZ-Two-Scales-RV / Gertler-Karadi-QE"
            }
        }
    }
}

/// A 2x4 sub-pixel canvas represented by Unicode Braille cells.
#[derive(Clone, Debug)]
pub struct BrailleCanvas {
    pub cols: usize,
    pub rows: usize,
    pub width: usize,
    pub height: usize,
    pixels: Vec<bool>,
    colors: Vec<[u8; 3]>,
}
impl BrailleCanvas {
    pub fn new(cols: usize, rows: usize) -> Self {
        let width = cols.saturating_mul(2);
        let height = rows.saturating_mul(4);
        Self {
            cols,
            rows,
            width,
            height,
            pixels: vec![false; width * height],
            colors: vec![[0; 3]; width * height],
        }
    }
    pub fn clear(&mut self) {
        self.pixels.fill(false);
        self.colors.fill([0; 3]);
    }
    pub fn set_pixel(&mut self, x: usize, y: usize, color: [u8; 3]) {
        if x < self.width && y < self.height {
            let i = y * self.width + x;
            self.pixels[i] = true;
            self.colors[i] = color;
        }
    }
    pub fn draw_line(
        &mut self,
        mut x0: isize,
        mut y0: isize,
        x1: isize,
        y1: isize,
        color: [u8; 3],
    ) {
        let dx = (x1 - x0).abs();
        let sx = if x0 < x1 { 1 } else { -1 };
        let dy = -(y1 - y0).abs();
        let sy = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;
        loop {
            if x0 >= 0 && y0 >= 0 {
                self.set_pixel(x0 as usize, y0 as usize, color);
            }
            if x0 == x1 && y0 == y1 {
                break;
            }
            let twice = 2 * err;
            if twice >= dy {
                err += dy;
                x0 += sx;
            }
            if twice <= dx {
                err += dx;
                y0 += sy;
            }
        }
    }
    pub fn plot_series(&mut self, xs: &[f64], ys: &[f64], y_min: f64, y_max: f64, color: [u8; 3]) {
        if xs.is_empty()
            || ys.is_empty()
            || xs.len() != ys.len()
            || self.width == 0
            || self.height == 0
        {
            return;
        }
        let x_min = xs.iter().copied().fold(f64::INFINITY, f64::min);
        let x_max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let dx = (x_max - x_min).abs().max(f64::EPSILON);
        let dy = (y_max - y_min).abs().max(f64::EPSILON);
        let canvas_width = self.width;
        let canvas_height = self.height;
        let point = |i: usize| {
            let x =
                (((xs[i] - x_min) / dx) * (canvas_width.saturating_sub(1) as f64)).round() as isize;
            let y_norm = ((ys[i] - y_min) / dy).clamp(0.0, 1.0);
            let y = ((1.0 - y_norm) * canvas_height.saturating_sub(1) as f64).round() as isize;
            (x, y)
        };
        let mut last = point(0);
        self.set_pixel(last.0.max(0) as usize, last.1.max(0) as usize, color);
        for i in 1..xs.len() {
            let next = point(i);
            self.draw_line(last.0, last.1, next.0, next.1, color);
            last = next;
        }
    }
    pub fn render_lines(&self, ansi_color: bool) -> Vec<String> {
        let mut lines = Vec::with_capacity(self.rows);
        for row in 0..self.rows {
            let mut line = String::new();
            for col in 0..self.cols {
                let mut mask = 0u32;
                let mut color = [0u8; 3];
                let mut active = false;
                for (py, dots) in BRAILLE_DOTS.iter().enumerate() {
                    for (px, dot) in dots.iter().enumerate() {
                        let x = col * 2 + px;
                        let y = row * 4 + py;
                        let idx = y * self.width + x;
                        if self.pixels[idx] {
                            mask |= *dot;
                            color = self.colors[idx];
                            active = true;
                        }
                    }
                }
                let ch = char::from_u32(0x2800 | mask).unwrap_or('\u{2800}');
                if ansi_color && active {
                    let _ = write!(
                        line,
                        "\x1b[38;2;{};{};{}m{}\x1b[0m",
                        color[0], color[1], color[2], ch
                    );
                } else {
                    line.push(ch);
                }
            }
            lines.push(line);
        }
        lines
    }
}

pub fn render_sparkline(values: &[f64], width: usize) -> String {
    if values.is_empty() || width == 0 {
        return String::new();
    }
    const LEVELS: [char; 8] = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
    let lo = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::INFINITY, f64::min);
    let hi = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    if !lo.is_finite() || !hi.is_finite() {
        return LEVELS[0].to_string().repeat(width);
    }
    let span = (hi - lo).abs().max(f64::EPSILON);
    (0..width)
        .map(|i| {
            let j = i * values.len() / width;
            let v = values[j.min(values.len() - 1)];
            let v = if v.is_finite() { v } else { lo };
            let level = (((v - lo) / span) * 7.0).round().clamp(0.0, 7.0) as usize;
            LEVELS[level]
        })
        .collect()
}

pub fn render_diverging_bar(
    value: f64,
    max_abs: f64,
    half_width: usize,
    ansi_color: bool,
) -> String {
    let n = if max_abs > 0.0 {
        ((value.abs() / max_abs) * half_width as f64).round() as usize
    } else {
        0
    };
    let n = n.min(half_width);
    let mut out = String::new();
    let color = if value < 0.0 {
        [255, 90, 90]
    } else {
        [0, 230, 140]
    };
    if ansi_color {
        let _ = write!(out, "\x1b[38;2;{};{};{}m", color[0], color[1], color[2]);
    }
    if value < 0.0 {
        out.push_str(&" ".repeat(half_width - n));
        out.push_str(&"━".repeat(n));
    } else {
        out.push_str(&" ".repeat(half_width));
    }
    out.push('│');
    if value >= 0.0 {
        out.push_str(&"━".repeat(n));
        out.push_str(&" ".repeat(half_width - n));
    } else {
        out.push_str(&" ".repeat(half_width));
    }
    if ansi_color {
        out.push_str("\x1b[0m");
    }
    out
}

fn visible_width(s: &str) -> usize {
    let mut count = 0;
    let mut escape = false;
    for ch in s.chars() {
        if escape {
            if ch == 'm' {
                escape = false;
            }
        } else if ch == '\x1b' {
            escape = true;
        } else {
            count += 1;
        }
    }
    count
}

/// Pad a dashboard cell by terminal-visible characters, ignoring ANSI color
/// sequences.  This keeps the 2x2 layout aligned in both plain and color mode.
fn pad_visible(s: &str, width: usize) -> String {
    let visible = visible_width(s);
    if visible >= width {
        return s.to_string();
    }
    let mut out = s.to_string();
    out.push_str(&" ".repeat(width - visible));
    out
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ScenarioKind {
    Normal,
    RatesPlus50Bp,
    EquityCrashMinus15Pct,
    CreditWideningPlus100Bp,
    StagflationPlus150Bp,
}
impl ScenarioKind {
    #[allow(non_upper_case_globals)]
    pub const RatesParallelPlus50Bp: Self = Self::RatesPlus50Bp;
    #[allow(non_upper_case_globals)]
    pub const VolShockPlus5Pts: Self = Self::EquityCrashMinus15Pct;
    #[allow(non_upper_case_globals)]
    pub const CreditSpreadWidening100Bp: Self = Self::CreditWideningPlus100Bp;
    pub fn next(self) -> Self {
        match self {
            Self::Normal => Self::RatesPlus50Bp,
            Self::RatesPlus50Bp => Self::EquityCrashMinus15Pct,
            Self::EquityCrashMinus15Pct => Self::CreditWideningPlus100Bp,
            Self::CreditWideningPlus100Bp => Self::StagflationPlus150Bp,
            Self::StagflationPlus150Bp => Self::Normal,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::Normal => "NORMAL",
            Self::RatesPlus50Bp => "RATES +50BP",
            Self::EquityCrashMinus15Pct => "EQUITY -15%",
            Self::CreditWideningPlus100Bp => "CREDIT +100BP",
            Self::StagflationPlus150Bp => "STAGFLATION +150BP",
        }
    }
}
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ChartMode {
    YieldAndInflationCurves,
    RoughHestonAndFxSmiles,
    XvaExposureProfile,
}
impl ChartMode {
    pub fn next(self) -> Self {
        match self {
            Self::YieldAndInflationCurves => Self::RoughHestonAndFxSmiles,
            Self::RoughHestonAndFxSmiles => Self::XvaExposureProfile,
            Self::XvaExposureProfile => Self::YieldAndInflationCurves,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::YieldAndInflationCurves => "Yield & Inflation Curves",
            Self::RoughHestonAndFxSmiles => "Rough Heston & FX Smiles",
            Self::XvaExposureProfile => "XVA Exposure Profile",
        }
    }
}

#[derive(Clone, Debug)]
pub struct WorkstationFrameData {
    pub tick: usize,
    pub scenario: ScenarioKind,
    pub chart_mode: ChartMode,
    pub eval_micros: f64,
    pub dag_recomputed_trades: usize,
    pub dag_total_trades: usize,
    pub dag_curves_rebuilt: bool,
    pub total_npv: f64,
    pub intraday_pnl: f64,
    pub eq_delta: f64,
    pub eq_vega: f64,
    pub fx_spot_delta: f64,
    pub inflation_par_delta_sum: f64,
    pub cva: f64,
    pub dva: f64,
    pub fva: f64,
    pub mva: f64,
    pub simm_initial_margin: f64,
    pub ir_dv01_buckets: Vec<(String, f64)>,
    pub credit_cs01_buckets: Vec<(String, f64)>,
    pub jump_to_default: f64,
    pub bermudan_sabr_metric: f64,
    pub sacc_swing_qp_metric: f64,
    pub slv_rough_worst_of_metric: f64,
    pub frtb_cvar_execution_metric: f64,
    pub vix_parisian_robust_metric: f64,
    pub double_heston_tarf_coco_hrp_metric: f64,
    pub jy_tba_bergomi_nco_metric: f64,
    pub cheyette_cmo_qrh_evar_metric: f64,
    pub asw_mf1f_ah_bess_ep_metric: f64,
    pub mtmxccs_colva_zabr_weather_enb_metric: f64,
    pub burnout_cmcds_quintic_afv_cs_schur_metric: f64,
    pub sinhou_mb_rbergomi_choi_afs_expectile_metric: f64,
    pub monotone_convex_clark_saccr_rvol_cartea_paleologo_metric: f64,
    pub krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric: f64,
    pub rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric: f64,
    pub mbs_rmv_mo_sig_garch_rie_gp_ierp_metric: f64,
    pub kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric: f64,
    pub mf_rates_wildcard_ctd_pivot_tarf_bergomi_metric: f64,
    pub cochrane_ftpl_g2_twist_essvi_offshore_cb_metric: f64,
    pub bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_metric: f64,
    pub veronesi_conv_lando_gen_gatheral_corridor_metric: f64,
    pub romer_taylor_damodaran_crp_christoffersen_var_metric: f64,
    pub fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric: f64,
    pub cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_metric: f64,
    pub rebonato_cms_leland_toft_ascot_pool_roll_rie_metric: f64,
    pub tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_metric: f64,
    pub hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_metric: f64,
    pub brigo_csa_ctd_callable_oas_jlt_corridor_lsm_metric: f64,
    pub brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_metric: f64,
    pub hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_metric: f64,
    pub brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_metric: f64,
    pub brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_metric: f64,
    pub andersen_cms_special_repo_cds_upfront_ssr_kva_crack_metric: f64,
    pub darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_metric: f64,
    pub henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_metric: f64,
    pub brigo_diff_swap_callable_muni_bkm_skew_break_cva_metric: f64,
    pub andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_metric: f64,
    pub litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_metric: f64,
    pub piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_metric: f64,
    pub gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_metric: f64,
    pub kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_metric: f64,
    pub jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_metric: f64,
    pub lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_metric: f64,
    pub jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_metric:
        f64,
    pub bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_metric:
        f64,
    pub andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_metric:
        f64,
    pub piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_metric:
        f64,
    pub mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_metric:
        f64,
    pub ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_metric:
        f64,
    pub smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_metric:
        f64,
    pub mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_metric:
        f64,
    pub hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_metric:
        f64,
    pub mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_metric:
        f64,
    pub filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_metric:
        f64,
    pub lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_metric:
        f64,
    pub cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_metric:
        f64,
    pub bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_metric:
        f64,
    pub hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_metric:
        f64,
    pub klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_metric:
        f64,
    pub blotter_rows: Vec<(String, f64, String)>,
}

fn range2(a: &[f64], b: &[f64]) -> (f64, f64) {
    let lo = a.iter().chain(b).copied().fold(f64::INFINITY, f64::min);
    let hi = a.iter().chain(b).copied().fold(f64::NEG_INFINITY, f64::max);
    let pad = ((hi - lo).abs() * 0.08).max(1e-6);
    (lo - pad, hi + pad)
}

fn range3(a: &[f64], b: &[f64], c: &[f64]) -> (f64, f64) {
    let lo = a
        .iter()
        .chain(b)
        .chain(c)
        .copied()
        .fold(f64::INFINITY, f64::min);
    let hi = a
        .iter()
        .chain(b)
        .chain(c)
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let pad = ((hi - lo).abs() * 0.08).max(1e-6);
    (lo - pad, hi + pad)
}

fn monotone_convex_clark_saccr_rvol_cartea_live_metric(
    m: &MarketSnapshot,
    credit_bump: f64,
) -> f64 {
    let c = FlatForwardCurve::new(m.ref_date, m.ois_quotes[0]);
    let hw = ql_curves::bootstrap_hagan_west_monotone_convex_curve(
        &ql_curves::HaganWestMonotoneConvexBootstrapSpec {
            reference_date: m.ref_date,
            pillar_times: vec![1.0, 2.0, 5.0],
            par_swap_quotes: m.ois_quotes.to_vec(),
            payment_frequency_per_year: 2,
            enforce_positivity: true,
            ameliorate_slopes: true,
            query_time: 2.0,
        },
    );
    let sw = ql_rates::price_piterbarg_quasi_gaussian_amortizing_swaption(
        &ql_rates::PiterbargQuasiGaussianAmortizingSwaptionSpec {
            option_expiry: 1.0,
            second_exercise_time: None,
            swap_times: vec![1.0, 2.0, 3.0],
            period_notionals: vec![1.0, 0.9],
            strike: m.ois_quotes[0],
            mean_reversion: 0.05,
            sigma0: 0.01 + m.eq_vol * 0.02,
            local_vol_skew_b: 0.0,
            cms_pay_delay: 0.0,
            is_payer: true,
        },
        &c,
        &c,
    );
    let fx = ql_fx::calibrate_clark_fx_delta_smile(&ql_fx::ClarkFxDeltaConventionSpec {
        spot: 1.10,
        r_dom: m.ois_quotes[0],
        r_for: m.fwd_quotes[0],
        expiry: 1.0,
        atm_vol: 0.12 + m.eq_vol * 0.1,
        rr_25d: 0.0,
        market_strangle_25d: 0.01,
        is_premium_included_delta: false,
        is_forward_delta: false,
        query_strike: 1.10,
    });
    let rv =
        ql_vol::evaluate_gatheral_gyongy_lsv_leverage(&ql_vol::GatheralLsvGyongyLeverageSpec {
            forward: m.eq_spot,
            strike: m.eq_spot,
            expiry: 1.0,
            svi_a: (m.eq_vol * m.eq_vol).max(1e-4),
            svi_b: 0.2,
            svi_rho: -0.3,
            svi_m: 0.0,
            svi_sigma: 0.2,
            dw_dt_atm: 0.04,
            heston_v0: 0.04,
            heston_theta: 0.04,
            heston_kappa: 1.0,
            heston_eta: 0.5,
            heston_rho: -0.7,
        });
    let sac = ql_xva::compute_green_saccr_lifetime_kva_mva(&ql_xva::GreenSaCcrLifetimeKvaMvaSpec {
        time_grid: vec![0.0, 1.0, 2.0],
        trade_notional: 1.0e6,
        trade_maturity: 5.0,
        initial_mtm: 0.0,
        mtm_volatility: 0.1,
        collateral_threshold: 0.0,
        initial_margin_ratio: 0.1,
        counterparty_pd: 0.02 + credit_bump,
        counterparty_lgd: 0.6,
        cva_risk_weight: 0.5,
        capital_ratio_requirement: 0.08,
        hurdle_rate_gamma_k: 0.1,
        funding_spread_sf: 0.01,
        mva_spread_sm: 0.01,
        risk_free_rate: m.ois_quotes[0],
    });
    let greg = ql_xva::allocate_gregory_netting_set_xva(&ql_xva::GregoryNettingSetAllocationSpec {
        trade_weights: vec![1.0, 1.0],
        trade_expected_mtms: vec![1.0, -1.0],
        trade_exposure_vols: vec![0.2, 0.2],
        correlation_matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        collateral_fraction: 0.0,
        horizon: 1.0,
        num_steps: 4,
        risk_free_rate: m.ois_quotes[0],
        cpty_hazard_rate: 0.02 + credit_bump,
        cpty_lgd: 0.6,
        own_hazard_rate: 0.01,
        own_lgd: 0.6,
        funding_spread: 0.01,
        colva_spread: 0.01,
    });
    let rough = ql_models::evaluate_bayer_friz_gatheral_markovian_rough(
        &ql_models::BayerFrizGatheralMarkovianRoughSpec {
            spot: m.eq_spot,
            strike: m.eq_spot,
            expiry: 1.0,
            rate: m.ois_quotes[0],
            div_yield: 0.0,
            xi0: (m.eq_vol * m.eq_vol).max(1e-4),
            eta: 1.0,
            rho: -0.7,
            hurst: 0.1,
            num_markov_factors: 4,
            is_call: true,
        },
    );
    let cli = ql_exotics::price_bouzoubaa_reverse_and_global_cliquet(
        &ql_exotics::BouzoubaaReverseAndGlobalCliquetSpec {
            notional: 1.0,
            rate: m.ois_quotes[0],
            div_yield: 0.0,
            reset_times: vec![0.5, 1.0],
            forward_vols: vec![m.eq_vol, m.eq_vol],
            forward_skews: vec![0.0, 0.0],
            initial_coupon_cmax: 0.0,
            local_floor: -0.1,
            local_cap: 0.1,
            global_floor: -0.2,
            global_cap: 0.2,
            is_reverse_cliquet: false,
        },
    );
    let plant = ql_commodities::price_burger_multifuel_tolling_power_plant(
        &ql_commodities::BurgerMultiFuelTollingPlantSpec {
            power_forwards: vec![80.0, 80.0],
            gas_forwards: vec![20.0, 20.0],
            coal_forwards: vec![15.0, 15.0],
            carbon_forwards: vec![50.0, 50.0],
            dispatch_times: vec![0.25, 0.75],
            capacity_mw: 100.0,
            period_hours: 1.0,
            gas_heat_rate: 2.0,
            gas_emission_factor: 0.4,
            gas_vom_cost: 2.0,
            coal_heat_rate: 2.5,
            coal_emission_factor: 0.8,
            coal_vom_cost: 1.0,
            vol_power: 0.2,
            vol_gas: 0.2,
            vol_coal: 0.2,
            vol_carbon: 0.2,
            rho_power_gas: 0.2,
            rho_power_coal: 0.2,
            rho_power_carbon: 0.2,
            rho_gas_coal: 0.2,
            start_up_cost: 100.0,
            min_load_fraction: 0.5,
            part_load_heat_rate_penalty: 0.1,
            risk_free_rate: m.ois_quotes[0],
        },
    );
    let exec = ql_portfolio::solve_cartea_jaimungal_alpha_execution(
        &ql_portfolio::CarteaJaimungalAlphaExecutionSpec {
            initial_inventory: 1.0,
            horizon: 1.0,
            num_steps: 8,
            temporary_impact_k: 1.0,
            permanent_impact_b: 0.1,
            inventory_penalty_phi: 0.2,
            terminal_penalty_a: 0.0,
            initial_alpha_signal: 0.1,
            alpha_mean_reversion_kappa: 1.0,
            alpha_volatility: 0.2,
        },
    );
    let paleo = ql_portfolio::solve_paleologo_factor_neutral_portfolio(
        &ql_portfolio::PaleologoFactorNeutralAlphaSpec {
            raw_alphas: vec![1.0, 2.0],
            factor_loadings_b: vec![vec![1.0], vec![1.0]],
            factor_covariance_omega: vec![vec![1.0]],
            sample_idio_variances: vec![0.1, 0.2],
            idio_shrinkage_xi: 0.5,
            idio_risk_aversion: 1.0,
            factor_risk_aversion: 1.0,
            factor_crowding_scores: vec![0.0, 0.0],
            crowding_penalty_lambda: 0.0,
        },
    );
    let aad = ql_dml::evaluate_glasserman_malliavin_checkpointed_aad(
        &ql_dml::GlassermanMalliavinCheckpointedMcSpec {
            spot: m.eq_spot,
            strike: m.eq_spot,
            rate: m.ois_quotes[0],
            div_yield: 0.0,
            base_vol: m.eq_vol,
            local_vol_skew: 0.0,
            expiry: 1.0,
            num_steps: 8,
            num_checkpoints: 2,
        },
    );
    hw.query_discount_factor.val()
        + sw.european_amortizing_npv.val()
        + fx.query_call_price.val()
        + rv.implied_vol.val()
        + sac.lifetime_kva.val()
        + greg.netted_total_xva.val()
        + rough.option_price.val()
        + cli.npv.val()
        + plant.tolling_plant_npv.val()
        + exec.expected_execution_cost_net_of_alpha.val()
        + paleo.expected_active_return.val()
        + aad.call_price.val()
}

pub struct DeskWorkstation {
    pub dag: PortfolioDag,
    pub base_market: MarketSnapshot,
    pub scenario: ScenarioKind,
    pub chart_mode: ChartMode,
    pub npv_history: Vec<f64>,
    trades: Vec<PortfolioTrade>,
    last_summary: PortfolioRiskSummary,
    tick: usize,
}
impl DeskWorkstation {
    pub fn new() -> Self {
        let ref_date = Date::from_ymd(2024, 1, 1).expect("valid workstation date");
        let market = MarketSnapshot {
            ref_date,
            ois_quotes: [0.02, 0.021, 0.022],
            fwd_quotes: [0.025, 0.026, 0.027],
            inflation_quotes: vec![("CPI".to_string(), 0.023)],
            heston_params: HestonParams {
                v0: 0.04,
                kappa: 1.5,
                theta: 0.04,
                sigma_v: 0.3,
                rho: -0.7,
            },
            eq_spot: 100.0,
            eq_vol: 0.2,
        };
        let program =
            ql_scripting::BytecodeProgram::compile(&ql_scripting::ScriptedProduct::autocallable(
                vec![0.25, 0.5, 0.75, 1.0],
                105.0,
                95.0,
                0.04,
                100.0,
                85.0,
                2.0,
            ));
        let trades = vec![
            PortfolioTrade::IrsSwap {
                id: 0,
                notional: 1_000_000.0,
                fixed_rate: 0.02,
                maturity_years: 5,
            },
            PortfolioTrade::IrsSwap {
                id: 1,
                notional: 1_000_000.0,
                fixed_rate: 0.021,
                maturity_years: 10,
            },
            PortfolioTrade::IrsSwap {
                id: 2,
                notional: 1_000_000.0,
                fixed_rate: 0.022,
                maturity_years: 20,
            },
            PortfolioTrade::HestonEuropeanOption {
                id: 3,
                strike: 100.0,
                expiry: 1.0,
                is_call: true,
                notional: 1.0,
            },
            PortfolioTrade::HestonEuropeanOption {
                id: 4,
                strike: 100.0,
                expiry: 2.0,
                is_call: false,
                notional: 1.0,
            },
            PortfolioTrade::ScriptedNote {
                id: 5,
                program,
                expiry: 1.0,
            },
        ];
        let mut dag = PortfolioDag::new();
        let summary = dag.evaluate_portfolio(&market, &trades);
        Self {
            dag,
            base_market: market,
            scenario: ScenarioKind::Normal,
            chart_mode: ChartMode::YieldAndInflationCurves,
            npv_history: vec![summary.total_npv],
            trades,
            last_summary: summary,
            tick: 0,
        }
    }
    pub fn set_scenario(&mut self, scenario: ScenarioKind) {
        self.scenario = scenario;
    }
    pub fn set_chart_mode(&mut self, mode: ChartMode) {
        self.chart_mode = mode;
    }
    fn scenario_market(&self) -> MarketSnapshot {
        let mut s = self.base_market.clone();
        match self.scenario {
            ScenarioKind::Normal => {
                s.eq_spot += self.tick as f64 * 0.0001;
            }
            ScenarioKind::RatesPlus50Bp => {
                s.ois_quotes.iter_mut().for_each(|x| *x += 0.005);
                s.fwd_quotes.iter_mut().for_each(|x| *x += 0.005);
            }
            ScenarioKind::EquityCrashMinus15Pct => {
                s.eq_spot *= 0.85;
                s.heston_params.v0 *= 1.5;
            }
            ScenarioKind::CreditWideningPlus100Bp => {}
            ScenarioKind::StagflationPlus150Bp => {
                s.ois_quotes.iter_mut().for_each(|x| *x += 0.0025);
                for (_, z) in &mut s.inflation_quotes {
                    *z += 0.015;
                }
            }
        }
        s
    }
    fn credit_metrics(&self, scenario: ScenarioKind) -> (Vec<(String, f64)>, f64, f64, f64) {
        let d = self.base_market.ref_date;
        let mats = [1.0, 2.0, 3.0, 5.0, 10.0];
        let quotes: Vec<_> = mats
            .iter()
            .enumerate()
            .map(|(i, &maturity_years)| CdsQuote {
                maturity_years,
                par_spread: 0.01
                    + i as f64 * 0.001
                    + if scenario == ScenarioKind::CreditWideningPlus100Bp {
                        0.01
                    } else {
                        0.0
                    },
            })
            .collect();
        let discount = DiscountCurve::new(
            d,
            mats.to_vec(),
            mats.iter().map(|t| (-0.02 * t).exp()).collect(),
        );
        let calibrated = CalibratedCreditCurve::bootstrap(d, discount.clone(), 0.4, &quotes)
            .expect("credit calibration");
        let cds = CreditDefaultSwap {
            side: ProtectionSide::Buyer,
            notional: 1_000_000.0,
            running_spread: 0.005,
            recovery_rate: 0.4,
            payment_times: (1..=20).map(|i| i as f64 / 4.0).collect(),
            accrual_fractions: vec![0.25; 20],
        };
        let (_, cs01) = calibrated.value_and_par_cs01(|curve| {
            let dc = FlatForwardCurve::new(d, ql_ad::TapeVar::from_f64(0.02));
            IsdaCdsEngine::price(&cds, &dc, curve).npv
        });
        let labels = ["1Y", "2Y", "3Y", "5Y", "10Y"];
        let buckets = cs01
            .iter()
            .zip(labels)
            .map(|(x, label)| (label.to_string(), *x))
            .collect();
        let jtd = CalibratedCreditCurve::jump_to_default(&cds, 0.0);
        let cds_npv = IsdaCdsEngine::price(&cds, &discount, &calibrated.curve).npv;
        (
            buckets,
            jtd,
            calibrated.curve.hazard_rates.iter().sum(),
            cds_npv,
        )
    }
    fn inflation_metrics(&self, scenario: ScenarioKind) -> (f64, f64) {
        let d = self.base_market.ref_date;
        let discount = DiscountCurve::new(
            d,
            vec![1.0, 2.0, 5.0],
            vec![(-0.02f64).exp(), (-0.04f64).exp(), (-0.1f64).exp()],
        );
        let bump = if scenario == ScenarioKind::StagflationPlus150Bp {
            0.015
        } else {
            0.0
        };
        let quotes = [(1.0, 0.02 + bump), (2.0, 0.022 + bump), (5.0, 0.025 + bump)];
        let mut f = [1.0; 12];
        for (i, factor) in f.iter_mut().enumerate().take(11) {
            *factor = 1.0 + 0.006 * (i as f64 - 5.5) / 5.5;
        }
        f[11] = 1.0 / f[..11].iter().product::<f64>();
        let seasonality = MultiplicativeSeasonality::new(f).expect("seasonality");
        let curve = CalibratedInflationCurve::bootstrap_with_seasonality(
            d,
            100.0,
            discount.clone(),
            &quotes,
            Some(seasonality),
        )
        .expect("inflation calibration");
        let swap = ZeroCouponInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.021,
            maturity_years: 3.0,
            target_month: Some(7),
            pay_fixed: true,
        };
        let (_, delta) = curve.value_and_par_inflation_delta(|c| {
            let dc = discount.to_ad_curve(
                &discount
                    .discounts
                    .iter()
                    .map(|&x| ql_ad::TapeVar::from_f64(x))
                    .collect::<Vec<_>>(),
            );
            swap.npv(&dc, c)
        });
        let sum: f64 = delta.iter().sum::<f64>() * 1e-4;
        let npv = swap.npv(&discount, &curve.curve).val();
        (sum, npv)
    }
    fn inflation_calibration(
        &self,
        scenario: ScenarioKind,
    ) -> (DiscountCurve<f64>, CalibratedInflationCurve) {
        let d = self.base_market.ref_date;
        let discount = DiscountCurve::new(
            d,
            vec![1.0, 2.0, 5.0, 10.0],
            vec![
                (-0.02f64).exp(),
                (-0.04f64).exp(),
                (-0.1f64).exp(),
                (-0.2f64).exp(),
            ],
        );
        let bump = if scenario == ScenarioKind::StagflationPlus150Bp {
            0.015
        } else {
            0.0
        };
        let quotes = [
            (1.0, 0.02 + bump),
            (2.0, 0.022 + bump),
            (5.0, 0.025 + bump),
            (10.0, 0.027 + bump),
        ];
        let mut factors = [1.0; 12];
        for (i, factor) in factors.iter_mut().enumerate().take(11) {
            *factor = 1.0 + 0.006 * (i as f64 - 5.5) / 5.5;
        }
        factors[11] = 1.0 / factors[..11].iter().product::<f64>();
        let seasonality = MultiplicativeSeasonality::new(factors).expect("seasonality");
        let curve = CalibratedInflationCurve::bootstrap_with_seasonality(
            d,
            100.0,
            discount.clone(),
            &quotes,
            Some(seasonality),
        )
        .expect("inflation calibration");
        (discount, curve)
    }
    fn active_curve_system(
        &self,
        market: &MarketSnapshot,
    ) -> Arc<ql_curves::CalibratedCurveSystem> {
        ql_portfolio::build_curve_system(market)
    }
    fn exposure_cube(&self) -> ql_xva::NpvCube<f64> {
        let d = self.base_market.ref_date;
        let market = self.scenario_market();
        let curve = FlatForwardCurve::new(d, market.ois_quotes[0]);
        let model = HullWhite1F {
            a: 0.05,
            sigma: 0.01 + 0.25 * market.eq_vol,
        };
        let swap = NettingSetSwap {
            notional: 1_000_000.0,
            fixed_rate: market.ois_quotes[1],
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accruals: vec![1.0; 5],
            is_payer: true,
        };
        let times: Vec<f64> = (1..=16).map(|k| k as f64 * 5.0 / 16.0).collect();
        simulate_hw1f_exposure_cube(&model, &curve, &[swap], &times, 64, 17)
    }
    fn xva_metrics(&self) -> (f64, f64, f64) {
        let d = self.base_market.ref_date;
        let market = self.scenario_market();
        let curve = FlatForwardCurve::new(d, market.ois_quotes[0]);
        let model = HullWhite1F {
            a: 0.05,
            sigma: 0.01,
        };
        let swap = NettingSetSwap {
            notional: 1_000_000.0,
            fixed_rate: market.ois_quotes[1],
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accruals: vec![1.0; 5],
            is_payer: true,
        };
        let cube =
            simulate_hw1f_exposure_cube(&model, &curve, &[swap], &[1.0, 2.0, 3.0, 4.0], 32, 17);
        let hazard = if matches!(self.scenario, ScenarioKind::CreditWideningPlus100Bp) {
            0.03
        } else {
            0.02
        };
        let report = compute_xva(&cube, &curve, None, hazard, 0.6, 0.015, 0.6, 0.005);
        (report.cva.val(), report.dva.val(), report.fva.val())
    }
    fn mva_metric(&self) -> f64 {
        let d = self.base_market.ref_date;
        let market = self.scenario_market();
        let curve = FlatForwardCurve::new(d, market.ois_quotes[0]);
        let model = HullWhite1F {
            a: 0.05,
            sigma: 0.01 + 0.25 * market.eq_vol,
        };
        let swap = NettingSetSwap {
            notional: 1_000_000.0,
            fixed_rate: market.ois_quotes[1],
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accruals: vec![1.0; 5],
            is_payer: true,
        };
        let times: Vec<f64> = (1..=12).map(|k| k as f64 * 5.0 / 12.0).collect();
        ql_xva::compute_chebyshev_dim_and_mva(&model, &curve, &[swap], &times, 64, 6, 0.01, 29).mva
    }
    fn krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_live_metric(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let r = m.ois_quotes[0] + credit_bump;
        let bg = ql_models::evaluate_bergomi_forward_start_and_var_strip(
            &ql_models::BergomiForwardStartAndVarStripSpec {
                params: ql_models::Bergomi2FParams {
                    xi0: m.eq_vol * m.eq_vol,
                    omega: 1.0,
                    theta: 0.5,
                    k1: 1.2,
                    k2: 2.4,
                    rho12: 0.2,
                    rho_s1: -0.5,
                    rho_s2: -0.3,
                },
                xi_long_run: m.eq_vol * m.eq_vol,
                xi_curve_decay: 1.0,
                forward_start_t1: 0.0,
                expiry_t2: 1.0,
                log_moneyness_k: 0.0,
            },
        )
        .forward_start_implied_vol
        .val();
        let levy = cos_levy_price(
            &VarianceGamma {
                sigma: m.eq_vol.max(1e-4),
                nu: 0.2,
                theta: -0.1,
            },
            m.eq_spot,
            m.eq_spot,
            r,
            0.01,
            1.0,
            ql_instruments::OptionType::Call,
            128,
        )
        .val();
        let cg = evaluate_cgmy_fractional_fft_and_pide(&ql_levy::CgmyFractionalFftAndPideSpec {
            spot: m.eq_spot,
            strikes: vec![m.eq_spot],
            risk_free_rate: r,
            dividend_yield: 0.01,
            expiry: 1.0,
            process: CgmyProcess {
                sigma: m.eq_vol.max(1e-4),
                c: 0.5,
                g: 5.0,
                m: 5.0,
                y: 0.7,
            },
            damping_alpha: 1.0,
            fourier_step_eta: 0.1,
            log_strike_step_lambda: 0.1,
            small_jump_cutoff_eps: 1e-4,
        })
        .frft_call_prices
        .first()
        .copied()
        .unwrap_or(0.0)
        .val();
        let mm = ql_portfolio::solve_gueant_lehalle_bounded_market_making(
            &ql_portfolio::GueantLehalleMarketMakingSpec {
                max_inventory_q: 2,
                current_inventory_q: 0,
                horizon: 1.0,
                num_time_steps: 4,
                volatility: m.eq_vol.max(1e-4),
                risk_aversion_gamma: 0.5,
                intensity_decay_k: 1.5,
                base_intensity_a: 2.0,
                terminal_inventory_penalty: 0.1,
                adverse_selection_zeta: 0.0,
            },
        )
        .current_optimal_bid_offset
        .val();
        let bl = ql_portfolio::solve_meucci_idzorek_black_litterman_and_plausibility(
            &ql_portfolio::MeucciIdzorekBlackLittermanSpec {
                covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
                market_weights: vec![0.6, 0.4],
                risk_aversion_lambda: 2.5,
                tau: 0.05,
                view_matrix_p: vec![vec![1.0, -1.0]],
                view_returns_q: vec![0.03],
                view_confidences: vec![0.65],
            },
        )
        .mahalanobis_plausibility_squared
        .val();
        let coco =
            price_despiegeleer_coco_and_mandatory_convertible(&DeSpiegeleerCoCoMandatoryConvSpec {
                notional: 100.0,
                spot: m.eq_spot,
                risk_free_rate: r,
                dividend_yield: 0.01,
                credit_spread: 0.02,
                maturity: 2.0,
                coupon_rate: 0.05,
                coupon_freq: 2,
                decs_lower_strike: 90.0,
                decs_upper_strike: 110.0,
                vol_lower_strike: m.eq_vol,
                vol_upper_strike: m.eq_vol,
                coco_cet1_trigger_spot: 70.0,
                coco_ponv_trigger_spot: 50.0,
                coco_conversion_floor_price: 60.0,
                coco_equity_vol: m.eq_vol,
            })
            .dual_trigger_coco_value
            .val();
        // Four remaining cross-domain pillars are evaluated on the same
        // live market snapshot (and credit bump) as the existing telemetry.
        let krd = analyze_tuckman_key_rate_durations_and_pca_hedge(&TuckmanKeyRateDurationSpec {
            nss_curve: NelsonSiegelSvensson {
                reference_date: m.ref_date,
                beta0: r + 0.01,
                beta1: -0.01,
                beta2: 0.02,
                beta3: -0.005,
                tau1: 1.5,
                tau2: 4.0,
            },
            cashflow_times: vec![1.0, 2.0, 3.0, 5.0],
            cashflow_amounts: vec![2.0, 2.0, 2.0, 102.0],
            key_rate_tenors: vec![1.0, 2.0, 5.0, 10.0],
            pca_decay_lambda: 0.25,
            benchmark_hedge_tenors: [2.0, 5.0, 10.0],
        })
        .dirty_price
        .val();
        let lpi =
            price_brigo_mercurio_yoy_and_lpi_market_model(&BrigoMercurioInflationMarketModelSpec {
                notional: 1_000_000.0,
                base_cpi: 100.0,
                payment_times: vec![1.0, 2.0, 3.0],
                forward_cpis: vec![102.0, 104.0, 106.0],
                nominal_forward_rates: vec![r; 3],
                nominal_discount_factors: vec![(-r).exp(), (-2.0 * r).exp(), (-3.0 * r).exp()],
                forward_cpi_vols: vec![m.eq_vol; 3],
                nominal_forward_vols: vec![m.eq_vol * 0.5; 3],
                inter_maturity_cpi_corr: 0.5,
                nominal_cpi_corr: 0.2,
                lpi_floor: -0.02,
                lpi_cap: 0.03,
                fixed_swap_rate: 0.02 + credit_bump,
            })
            .fimm_convexity_adjusted_lpi_swap_npv
            .val();
        let cdo =
            price_okane_bespoke_tranche_and_student_t_basket(&OKaneBespokeBaseCorrelationSpec {
                index_detachments: vec![0.03, 0.07, 0.10, 0.15],
                index_base_correlations: vec![0.35, 0.40, 0.45, 0.50],
                index_expected_loss: 0.04 + credit_bump,
                bespoke_hazard_rates: vec![0.02 + credit_bump; 5],
                bespoke_recovery_rate: 0.4,
                bespoke_attachment: 0.03,
                bespoke_detachment: 0.10,
                strike_mapping_exponent_gamma: 1.0,
                student_t_dof_nu: 6.0 + 10.0 * credit_bump,
                maturity: 5.0,
                risk_free_rate: r,
            })
            .student_t_copula_senior_tranche_loss
            .val();
        let fx = evaluate_wystup_fx_accumulator_and_dcd(&WystupFxStructuredProductSpec {
            spot: m.eq_spot,
            accumulator_strike: m.eq_spot * 0.95,
            knockout_barrier: m.eq_spot * 1.10,
            downside_gearing: 2.0,
            domestic_rate: r,
            foreign_rate: 0.01,
            volatility: m.eq_vol.max(1e-4),
            fixing_times: vec![0.25, 0.5, 0.75, 1.0],
            guaranteed_fixings: 1,
            dcd_conversion_strike: m.eq_spot,
            dcd_tenor: 1.0,
        })
        .total_accumulator_npv
        .val();
        let pde = ql_pde::solve_rannacher_crank_nicolson_sinh_pde(&ql_pde::RannacherSinhPdeSpec {
            spot: m.eq_spot,
            strike: m.eq_spot,
            volatility: m.eq_vol.max(1e-4),
            risk_free_rate: r,
            dividend_yield: 0.01,
            expiry: 1.0,
            is_call: true,
            n_space: 31,
            n_time: 8,
            rannacher_steps: 2,
            sinh_concentration: 0.15,
        })
        .pde_price
        .val();
        let mot =
            ql_ad::solve_mot_super_replication_bounds_and_ift(&ql_ad::MotSuperReplicationIftSpec {
                support_grid: vec![0.8, 1.0, 1.2],
                marginal_t1_probs: vec![0.2, 0.5, 0.3],
                marginal_t2_probs: vec![0.25, 0.5, 0.25],
                forward_start_strike: 0.95 + 0.001 * credit_bump,
                entropic_regularization_eps: 0.02,
            })
            .entropic_mot_price
            .val();
        let exotic = ql_dml::evaluate_ift_calibrated_exotic_quote_hedge(
            &ql_dml::IftCalibratedExoticHedgeSpec {
                spot: m.eq_spot,
                risk_free_rate: r,
                expiry: 1.0,
                calibration_strikes: [0.9 * m.eq_spot, m.eq_spot, 1.1 * m.eq_spot],
                market_implied_vols: [m.eq_vol * 1.1, m.eq_vol, m.eq_vol * 1.05],
                exotic_strike: 1.05 * m.eq_spot,
                exotic_power_p: 1.2,
            },
        )
        .exotic_option_price
        .val();
        let revolve = ql_ad::execute_savine_revolve_checkpointed_adjoint(
            &ql_ad::SavineRevolveCheckpointSpec {
                num_time_steps: 8,
                num_checkpoints: 2,
                initial_spot: m.eq_spot,
                drift_rate: r,
                volatility: m.eq_vol.max(1e-4),
                mean_reversion_kappa: 0.5,
                dt: 0.125,
                strike: m.eq_spot,
            },
        )
        .terminal_option_payoff_value
        .val();
        monotone_convex_clark_saccr_rvol_cartea_live_metric(m, credit_bump)
            + bg
            + levy
            + cg
            + mm
            + bl
            + coco
            + krd * 1e-4
            + lpi * 1e-6
            + cdo
            + fx
            + pde
            + mot
            + exotic
            + revolve
    }
    pub fn mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let ie_growth = (0.02 + 0.25 * (inf - 0.023)).clamp(0.005, 0.045);
        let mbs_cpr = (0.06 + 0.5 * (inf - 0.023)).clamp(0.02, 0.20);
        let mbs = ql_fixed_income::analyze_davidson_levin_mbs_io_po_strips(
            &ql_fixed_income::DavidsonLevinIoPoStripMbsSpec {
                notional: 100.0,
                wac_annual: 0.065,
                pass_through_coupon_annual: 0.06,
                current_mortgage_rate: m.ois_quotes[0] + credit_bump,
                discount_rate_annual: m.ois_quotes[0],
                remaining_months: 24,
                loan_age_months: 12,
                turnover_base_cpr: mbs_cpr,
                max_refi_cpr: 0.4,
                refi_slope: 80.0,
                refi_threshold_spread: 0.005,
                burnout_exponent: 1.0,
                rate_bump_bps: 1.0,
            },
        );
        let bond = ql_credit::analyze_duffie_singleton_rmv_rfv_and_cds_basis(
            &ql_credit::DuffieSingletonRecoveryAndCdsBasisSpec {
                face_value: 100.0,
                coupon_rate_annual: 0.05,
                payment_freq_per_year: 2,
                maturity_years: 3.0,
                risk_free_rate: m.ois_quotes[0],
                hazard_rate: 0.02 + credit_bump,
                recovery_rate: 0.4,
                market_bond_price: 98.0,
                repo_funding_spread: 0.005,
            },
        );
        let mo = ql_xva::compute_crepey_bielecki_marshall_olkin_bilateral_xva(
            &ql_xva::CrepeyBieleckiMarshallOlkinBilateralXvaSpec {
                times: vec![0.5, 1.0],
                expected_positive_exposure: vec![1.0, 1.0],
                expected_negative_exposure: vec![0.5, 0.5],
                risk_free_rate: m.ois_quotes[0],
                lambda_institution_idio: 0.01,
                lambda_counterparty_idio: 0.02,
                lambda_common_shock_12: 0.01 + credit_bump,
                recovery_institution: 0.4,
                recovery_counterparty: 0.4,
                recovery_counterparty_simultaneous: 0.3,
                gap_exposure_multiplier: 0.2,
                funding_borrowing_spread: 0.01,
                funding_lending_spread: 0.005,
            },
        );
        let cb = ql_fixed_income::price_sirbu_shreve_dilution_and_makewhole_convertible(
            &ql_fixed_income::SirbuShreveDilutionMakeWholeConvertibleSpec {
                spot: m.eq_spot,
                volatility: m.eq_vol,
                risk_free_rate: m.ois_quotes[0],
                dividend_yield: 0.0,
                credit_spread: 0.02 + credit_bump,
                maturity_years: 2.0,
                par_value: 100.0,
                coupon_rate_annual: 0.03,
                conversion_ratio: 0.8,
                existing_shares_millions: 100.0,
                convertible_issue_count_millions: 10.0,
                debt_relief_retention_fraction: 0.0,
                call_price: 105.0,
                soft_call_trigger_spot: 110.0,
                make_whole_max_bonus_shares: 0.2,
                make_whole_time_decay: 0.2,
                make_whole_moneyness_decay: 0.3,
                tree_steps: 8,
            },
        );
        let sig = ql_models::evaluate_bayer_oberhauser_signature_option_pricing(
            &ql_models::BayerHorvathOberhauserPathSignatureSpec {
                time_grid: vec![0.0, 0.5, 1.0],
                sample_paths: vec![vec![m.eq_spot, m.eq_spot * 1.01, m.eq_spot]],
                strike: m.eq_spot,
                discount_factor: (-m.ois_quotes[0]).exp(),
                ridge_lambda: 1e-6,
            },
        );
        let hn = ql_vol::evaluate_heston_nandi_and_tsay_dcc_garch(
            &ql_vol::HestonNandiTsayDccGarchSpec {
                spot: m.eq_spot,
                strike: m.eq_spot,
                daily_risk_free_rate: m.ois_quotes[0] / 252.0,
                days_to_maturity: 30,
                initial_daily_variance: m.eq_vol * m.eq_vol / 252.0,
                omega: 1e-5,
                alpha: 0.05,
                beta: 0.9,
                gamma_asymmetry: 0.2,
                asset1_returns: vec![0.01],
                asset2_returns: vec![0.008],
                dcc_a: 0.02,
                dcc_b: 0.95,
            },
        );
        let pot = ql_portfolio::evaluate_mcneil_embrechts_pot_gpd_and_christoffersen_backtest(
            &ql_portfolio::McNeilEmbrechtsPotGpdBacktestSpec {
                historical_losses: vec![1.0, 2.0, 3.0, 4.0, 5.0],
                threshold_u: 1.0,
                confidence_alpha: 0.95,
                spectral_risk_aversion_k: 2.0,
                backtest_var_threshold: 3.0,
            },
        );
        let rie = ql_portfolio::solve_potters_bouchaud_rie_nonlinear_shrinkage(
            &ql_portfolio::PottersBouchaudRieCovarianceSpec {
                sample_covariance: vec![vec![0.04 + credit_bump, 0.01], vec![0.01, 0.09]],
                num_observations_t: 100,
            },
        );
        let pin = ql_portfolio::evaluate_ohara_ekop_kyle_microstructure_and_etf_band(
            &ql_portfolio::OharaEkopKyleMicrostructureSpec {
                fundamental_value_std: 1.0,
                noise_trader_volume_std: 2.0,
                alpha_info_event_prob: 0.2,
                delta_bad_news_prob: 0.4,
                mu_informed_arrival_rate: 1.0,
                epsilon_uninformed_buy_rate: 2.0,
                epsilon_uninformed_sell_rate: 2.0,
                fundamental_low: 95.0,
                fundamental_high: 105.0,
                position_value: 100.0,
                return_vol_daily: 0.02,
                z_score_alpha: 2.0,
                mean_relative_spread: 0.01,
                std_relative_spread: 0.002,
                etf_price: m.eq_spot,
                etf_inav: 100.0,
                ap_creation_fee_bps: 5.0,
                basket_spread_bps: 2.0,
                etf_spread_bps: 3.0,
            },
        );
        let gp = ql_portfolio::solve_garleanu_pedersen_dynamic_aim_portfolio(
            &ql_portfolio::GarleanuPedersenDynamicAimPortfolioSpec {
                covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
                current_weights: vec![0.5, 0.5],
                slow_alpha_vector: vec![0.03, 0.02],
                fast_alpha_vector: vec![0.01, 0.02],
                slow_mean_reversion_phi: 0.1,
                fast_mean_reversion_phi: 0.5,
                risk_aversion_gamma: 2.0,
                transaction_cost_kappa: 0.1,
                discount_rate_r: m.ois_quotes[0],
            },
        );
        let br = ql_portfolio::evaluate_bacon_carino_brinson_and_karnosky_singer_attribution(
            &ql_portfolio::BaconCarinoMultiPeriodBrinsonSpec {
                portfolio_weights_by_period: vec![vec![0.5, 0.5]],
                benchmark_weights_by_period: vec![vec![0.5, 0.5]],
                portfolio_returns_by_period: vec![vec![0.1, 0.02]],
                benchmark_returns_by_period: vec![vec![0.08, 0.01]],
                fx_spot_returns_by_sector: vec![0.01],
                foreign_cash_rates_by_sector: vec![0.02],
                domestic_cash_rate: m.ois_quotes[0],
            },
        );
        let ie = ql_ad::solve_penman_damodaran_riv_aeg_and_implied_erp_ift(
            &ql_ad::PenmanDamodaranResidualIncomeIerpSpec {
                initial_book_value: 100.0,
                forecast_roe_schedule: vec![0.12, 0.11],
                payout_ratio_schedule: vec![0.4, 0.4],
                terminal_roe: 0.1,
                terminal_growth_rate: ie_growth,
                risk_free_rate: m.ois_quotes[0],
                equity_beta: 1.1,
                observed_market_price: 110.0,
            },
        );
        mbs.io_strip_price
            + bond.raw_cds_bond_basis
            + mo.total_bilateral_xva
            + cb.make_whole_enhanced_price
            + sig.signature_linear_functional_price
            + hn.heston_nandi_call_price
            + pot.pot_gpd_expected_shortfall
            + rie.rie_portfolio_variance
            + pin.ekop_pin
            + gp.next_period_optimal_weights[0]
            + br.carino_linked_selection_total
            + ie.implied_erp
    }
    fn rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_nco_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let rk = ql_curves::calibrate_filipovic_rkhs_smooth_yield_curve(
            &ql_curves::FilipovicRkhsYieldCurveSpec {
                cashflow_times: vec![1.0, 2.0],
                instrument_cashflow_matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
                instrument_dirty_prices: vec![m.ois_quotes[0], m.ois_quotes[1]],
                decay_alpha: 1.0,
                ridge_penalty_lambda: 1e-8,
                query_times: vec![1.0, 2.0],
            },
        );
        let fmm = ql_rates::price_lyashenko_mercurio_fmm_caplet_and_swaption(
            &ql_rates::LyashenkoMercurioFmmSpec {
                notional: 1.0,
                exercise_time_t0: 1.0,
                period_accruals: vec![0.5, 0.5],
                forward_rates: vec![m.ois_quotes[1], m.ois_quotes[1]],
                discount_factors: vec![0.98, 0.96],
                fmm_volatilities: vec![m.eq_vol, m.eq_vol],
                displacement_shift: 0.0,
                inter_period_correlation_decay: 1.0,
                strike: m.ois_quotes[1],
                is_payer_or_caplet: true,
            },
        );
        let vv = ql_fx::evaluate_castagna_full_vanna_volga_hedge_and_barrier(
            &ql_fx::CastagnaFullVannaVolgaHedgeSpec {
                spot: m.eq_spot,
                r_dom: m.ois_quotes[0],
                r_for: 0.01,
                expiry: 1.0,
                atm_vol: m.eq_vol,
                rr_25d: 0.0,
                bf_25d: 0.0,
                target_strike: m.eq_spot,
                knockout_barrier: m.eq_spot * 1.2,
                unwind_attenuation_factor: 0.5,
                is_up_and_out: true,
                is_call: true,
            },
        );
        let lsv = ql_vol::evaluate_guyon_markovian_projection_lsv_smile(
            &ql_vol::GuyonMarkovianProjectionLsvSpec {
                forward: 100.0,
                strike: 100.0,
                expiry: 1.0,
                forward_start_t1: 0.25,
                atm_vol: m.eq_vol,
                market_atm_skew: -0.1,
                market_atm_curvature: 0.1,
                sv_v0: m.eq_vol * m.eq_vol,
                sv_kappa: 1.0,
                sv_eta: 0.3,
                sv_rho: -0.5,
                sv_mixing_weight_eta: 0.5,
            },
        );
        let jlt = ql_credit::evaluate_jlt_rating_migration_and_crcva(
            &ql_credit::JarrowLandoTurnbullRatingMigrationSpec {
                initial_rating_index: 0,
                physical_generator_4x4: [
                    [-0.1, 0.08, 0.01, 0.01],
                    [0.02, -0.08, 0.04, 0.02],
                    [0.0, 0.02, -0.05, 0.03],
                    [0.0, 0.0, 0.0, 0.0],
                ],
                risk_premium_multipliers: [1.0, 1.1, 1.2],
                recovery_rate: 0.4,
                risk_free_rate: m.ois_quotes[0],
                maturity: 5.0,
                base_coupon_rate: 0.03,
                hy_downgrade_stepup_spread: 0.01,
                counterparty_epe_notional: 1.0,
            },
        );
        let ctd = ql_xva::compute_brigo_pallavicini_ctd_collateral_and_closeout_xva(
            &ql_xva::BrigoPallaviciniCtdCollateralXvaSpec {
                notional_mtm: 1.0,
                maturity: 1.0,
                num_steps: 8,
                base_ois_rate: m.ois_quotes[0],
                csa_spread_ccy1: 0.001,
                csa_spread_ccy2: 0.002,
                csa_spread_ccy3: 0.003,
                csa_mean_reversion: 1.0,
                csa_vol_ccy1: 0.01,
                csa_vol_ccy2: 0.01,
                csa_vol_ccy3: 0.01,
                csa_correlation: 0.2,
                cpty_hazard_rate: 0.02,
                cpty_lgd: 0.6,
                own_hazard_rate: 0.01,
                own_lgd: 0.6,
                exposure_volatility: 0.1,
            },
        );
        let qrh = ql_models::evaluate_gatheral_rosenbaum_joint_spx_vix_qrh(
            &ql_models::GatheralRosenbaumJointSpxVixQrhSpec {
                params: ql_models::QuadraticRoughHestonParams {
                    hurst: 0.1,
                    nu: 0.3,
                    lambda: 1.0,
                    a: 1.0,
                    b: 0.0,
                    c: 0.04,
                    z0: 0.0,
                },
                spot: 100.0,
                spx_strike: 100.0,
                vix_strike: 0.2,
                expiry: 1.0,
                rate: m.ois_quotes[0],
            },
        );
        let par = ql_exotics::price_chesney_yor_cumulative_parisian_and_rainbow(
            &ql_exotics::ChesneyYorCumulativeParisianRainbowSpec {
                spot_1: 100.0,
                spot_2: 102.0,
                strike: 100.0,
                par_barrier: 120.0,
                excursion_window_d: 0.1,
                vol_1: m.eq_vol,
                vol_2: m.eq_vol,
                div_1: 0.0,
                div_2: 0.0,
                correlation: 0.5,
                risk_free_rate: m.ois_quotes[0],
                expiry: 1.0,
                is_up_and_out: true,
            },
        );
        let aid = ql_commodities::price_aid_carmona_structural_electricity_stack(
            &ql_commodities::AidCarmonaStructuralElectricityStackSpec {
                gas_forward: m.fwd_quotes[0] * 1000.0,
                gas_vol: 0.3,
                base_heat_rate: 2.0,
                total_capacity_mw: 100.0,
                baseline_capacity_anchor_mw: 50.0,
                expected_load_mw: 80.0,
                expected_renewable_gen_mw: 20.0,
                residual_load_vol_mw: 10.0,
                stack_convexity_alpha: 1.0,
                negative_price_floor_depth: 5.0,
                curtailment_beta: 1.0,
                load_gas_correlation: 0.2,
                power_option_strike: 50.0,
                spike_threshold_price: 100.0,
                risk_free_rate: m.ois_quotes[0],
                expiry: 1.0,
            },
        );
        let imp = ql_portfolio::solve_bouchaud_webster_square_root_and_propagator_execution(
            &ql_portfolio::BouchaudWebsterPropagatorImpactSpec {
                metaorder_size_q: 1.0,
                daily_volume_vd: 10.0,
                daily_volatility_sigma: 0.2,
                y_prefactor: 1.0,
                propagator_exponent_beta: 0.5,
                propagator_scale_g0: 1.0,
                regularization_tau0: 0.01,
                temporary_impact_eta: 0.1,
                risk_aversion_lambda: 0.1,
                horizon: 1.0,
                num_slices: 8,
            },
        );
        let nco = ql_portfolio::solve_lopez_de_prado_denoised_nco_allocation(
            &ql_portfolio::LopezDePradoNcoDenoisedSpec {
                covariance: vec![
                    vec![0.04 + credit_bump, 0.01],
                    vec![0.01, 0.09 + credit_bump],
                ],
                expected_returns: vec![0.05, 0.06],
                sample_ratio_q: 2.0,
                denoising_alpha: 0.5,
                num_clusters: 2,
            },
        );
        let dml = ql_dml::evaluate_glasserman_savine_hybrid_pw_lrm_and_dml(
            &ql_dml::GlassermanSavineHybridPwLrmSpec {
                spot: 100.0,
                strike: 100.0,
                rate: m.ois_quotes[0],
                div_yield: 0.0,
                volatility: m.eq_vol,
                expiry: 1.0,
                digital_smoothing_eps: 0.5,
                dml_differential_weight_gamma: 1.0,
            },
        );
        rk.query_zero_rates[0].val()
            + fmm.backward_looking_rfr_caplet_npv.val()
            + vv.exact_vv_vanilla_price.val()
            + lsv.lsv_reconstructed_implied_vol.val()
            + jlt.rating_migration_adjusted_crcva.val()
            + ctd.ctd_collateral_switch_option_value.val()
            + qrh.vix_futures_fair_value.val()
            + par.rainbow_call_on_max.val()
            + aid.structural_power_forward.val()
            + imp.optimal_propagator_cost.val()
            + nco.nco_portfolio_volatility.val()
            + dml.hybrid_pw_lrm_digital_delta.val()
            + credit_bump
            + Self::krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_live_metric(m, credit_bump) * 1e-3
    }
    pub fn rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_nco_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let bump = if matches!(self.scenario, ScenarioKind::CreditWideningPlus100Bp) {
            0.01
        } else {
            0.0
        };
        Self::rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_nco_live_metric_from_market(&market, bump)
    }
    pub fn kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let r = m.ois_quotes[0];
        let spot = m.eq_spot;
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let kal = ql_curves::evaluate_tsay_kalman_rts_state_space_curve(
            &ql_curves::TsayKalmanRtsStateSpaceCurveSpec {
                transition_f: [[0.9, 0.1], [0.0, 0.95]],
                state_drift_c: [0.0, 0.0],
                process_noise_q: [[0.01, 0.0], [0.0, 0.01]],
                measurement_h: [[1.0, 0.0], [0.0, 1.0]],
                measurement_noise_r: [0.02, 0.02],
                initial_state_mean: [r, spot],
                initial_state_cov: [[1.0, 0.0], [0.0, 1.0]],
                observations: vec![[spot, r], [spot * 1.001, r + 0.001]],
            },
        );
        let vecm = ql_portfolio::solve_johansen_tsay_vecm_cointegration_basket(
            &ql_portfolio::JohansenTsayVecmCointegrationSpec {
                price_series: (0..8)
                    .map(|i| {
                        [
                            100.0 + i as f64,
                            100.0 + i as f64 + 0.5 * (1.0 - (0.9_f64).powi(i)),
                        ]
                    })
                    .collect(),
            },
        );
        let hj = ql_portfolio::evaluate_cochrane_hansen_jagannathan_sdf_and_campbell_shiller(
            &ql_portfolio::CochraneHansenJagannathanSdfSpec {
                mean_excess_returns: [0.03, 0.04],
                excess_return_cov: [[0.04, 0.01], [0.01, 0.09]],
                risk_free_rate: r,
                candidate_sdf_pricing_errors: [0.001, 0.002],
                steady_state_log_dp_ratio: inf.ln(),
                var1_companion_matrix: [[0.9, 0.1], [0.0, 0.8]],
                state_innovation_cov: [[0.02, 0.005], [0.005, 0.03]],
            },
        );
        let ifrs = ql_credit::compute_ifrs9_cecl_multistage_expected_credit_loss(
            &ql_credit::Ifrs9CeclMultiStageEclSpec {
                annual_ttc_pds: [
                    0.02 + credit_bump * 0.2,
                    0.04 + credit_bump * 0.2,
                    0.08 + credit_bump * 0.2,
                ],
                annual_eads: [100.0, 98.0, 95.0],
                base_lgd: 0.45,
                lgd_macro_sensitivity: 0.05,
                vasicek_asset_correlation: 0.2,
                effective_interest_rate: r,
                scenario_weights: [0.3, 0.5, 0.2],
                scenario_macro_shocks_z: [-1.0, 0.0, 1.0],
                origination_lifetime_pd: 0.04,
                sicr_pd_multiplier_threshold: 1.5,
            },
        );
        let cpt = ql_portfolio::evaluate_wakker_tversky_kahneman_cumulative_prospect_utility(
            &ql_portfolio::WakkerTverskyKahnemanCptSpec {
                scenario_returns: [-0.10, -0.02, 0.03, 0.12],
                scenario_probabilities: [0.1, 0.2, 0.4, 0.3],
                reference_return: 0.0,
                alpha_gain: 0.88,
                beta_loss: 0.88,
                loss_aversion_lambda: 2.25,
                prelec_gamma_gain: 0.7,
                prelec_gamma_loss: 0.8,
                prelec_delta: 1.0,
            },
        );
        let has = ql_portfolio::evaluate_hasbrouck_information_share_and_dssw_limits_to_arb(
            &ql_portfolio::HasbrouckDsswPriceDiscoverySpec {
                alpha_venue1: 0.3,
                alpha_venue2: 0.2,
                sigma_venue1: vol,
                sigma_venue2: vol * 1.1,
                venue_innovation_corr: 0.2,
                fundamental_dividend: 2.0,
                risk_free_rate: r,
                noise_trader_share_mu: 0.2,
                noise_trader_sentiment_mean: 0.01,
                noise_trader_sentiment_var: 0.02,
                arbitrageur_risk_aversion_gamma: 2.0,
            },
        );
        let fv = ql_rates::price_fong_vasicek_stochastic_vol_bond(
            &ql_rates::FongVasicekStochVolShortRateSpec {
                initial_short_rate_r0: r,
                initial_variance_v0: vol * vol,
                rate_mean_reversion_alpha: 0.5,
                long_run_rate_r_bar: r,
                vol_mean_reversion_gamma: 1.0,
                long_run_variance_v_bar: vol * vol,
                vol_of_vol_xi: 0.2,
                correlation_rho: -0.3,
                maturity_tau: 2.0,
            },
        );
        let lando = ql_credit::evaluate_lando_cox_affine_default_and_liquidity_decomposition(
            &ql_credit::LandoCoxAffineDefaultSpec {
                risk_free_rate: r,
                base_default_hazard: 0.01 + credit_bump,
                initial_credit_intensity_y0: 0.02,
                credit_kappa: 0.5,
                credit_theta: 0.02,
                credit_sigma: 0.1,
                initial_liquidity_intensity_z0: 0.01,
                liquidity_kappa: 0.6,
                liquidity_theta: 0.01,
                liquidity_sigma: 0.08,
                recovery_rate: 0.4,
                maturity: 5.0,
            },
        );
        let times: Vec<f64> = (0..8).map(|i| i as f64 + 1.0).collect();
        let logs: Vec<f64> = times
            .iter()
            .map(|&t| 4.0 - 0.2 * (10.0 - t).powf(0.5) + 0.01 * (8.0 * (10.0 - t).ln()).cos())
            .collect();
        let lppl = ql_portfolio::calibrate_sornette_lppl_bubble_and_crash_hazard(
            &ql_portfolio::SornetteLpplBubbleCrashSpec {
                observation_times: times,
                log_prices: logs,
                critical_time_tc: 10.0,
                power_law_exponent_m: 0.5,
                log_periodic_omega: 8.0,
                expected_crash_drop_kappa: 0.2,
                ridge_penalty: 1e-8,
            },
        );
        let cb = ql_fixed_income::analyze_thorp_kassouf_calamos_convertible_arbitrage(
            &ql_fixed_income::ThorpKassoufCalamosConvertibleArbSpec {
                stock_price: spot,
                conversion_ratio: 0.8,
                convertible_market_price: 95.0,
                face_value: 100.0,
                annual_coupon_rate: 0.03,
                stock_dividend_yield: 0.01,
                risk_free_rate: r,
                credit_spread: 0.03 + credit_bump,
                short_stock_rebate_rate: 0.01,
                implied_volatility: vol,
                realized_volatility: vol * 1.1,
                maturity: 3.0,
            },
        );
        let causal = ql_portfolio::evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(
            &ql_portfolio::LopezDePradoCausalFactorSpec {
                true_causal_beta_xy: 0.5,
                confounder_to_factor_axz: 0.8,
                confounder_to_return_byz: 0.6,
                factor_to_collider_ccx: 0.4,
                return_to_collider_ccy: 0.5,
                confounder_var: 1.0,
                factor_idiosyncratic_var: 0.5,
                return_idiosyncratic_var: 0.7,
                collider_idiosyncratic_var: 0.4,
            },
        );
        let mal = ql_dml::evaluate_fournie_lasry_lions_malliavin_bel_greeks(
            &ql_dml::FournieLasryLionsMalliavinBelSpec {
                spot,
                strike: spot,
                rate: r,
                div_yield: 0.0,
                vol,
                expiry: 1.0,
                localization_bandwidth: 0.5,
            },
        );
        kal.innovation_log_likelihood.val()
            + vecm.johansen_trace_stat_r0.val() * 1e-3
            + hj.discount_rate_news_variance_share.val()
            + ifrs.blended_ifrs9_ecl_provision.val() * 1e-2
            + cpt.cumulative_prospect_utility.val()
            + has.efficient_price_innovation_var.val()
            + fv.zero_coupon_bond_price.val()
            + lando.risky_corporate_bond_price.val()
            + lppl.terminal_crash_hazard_rate.val()
            + cb.theoretical_convertible_value.val() * 1e-2
            + causal.kitchen_sink_collider_conditioned_beta.val()
            + mal.call_malliavin_delta.val()
            + inf * 10.0
            + credit_bump
    }
    pub fn cochrane_ftpl_g2_twist_essvi_offshore_cb_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let r = m.ois_quotes[0] + credit_bump;
        let spot = m.eq_spot;
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let ftpl = ql_inflation::evaluate_cochrane_ftpl_and_maturity_inflation_smoothing(
            &ql_inflation::CochraneFiscalTheoryPriceLevelSpec {
                debt_maturities_years: vec![1.0, 5.0, 10.0],
                nominal_debt_face_values: vec![30.0, 30.0, 40.0],
                nominal_zero_yields: vec![r, r + 0.002, r + 0.004],
                real_zero_yields: vec![r - inf; 3],
                primary_real_surpluses: vec![5.0, 5.2, 5.4],
                terminal_surplus_growth_rate: 0.02,
                fiscal_surplus_shock_fraction: -0.05,
                expected_future_inflation_rise: inf,
                inflation_risk_beta: 0.4,
                inflation_volatility: vol,
                tips_liquidity_spread: 0.002,
            },
        );
        let g2 = ql_rates::price_brigo_mercurio_g2_twist_and_swaption(
            &ql_rates::BrigoMercurioG2TwistSwaptionSpec {
                a: 0.1,
                sigma: vol,
                b: 0.3,
                eta: vol * 0.8,
                rho: -0.4,
                flat_zero_rate: r,
                option_expiry: 1.0,
                swap_payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
                fixed_strike: r,
                notional: 1.0,
            },
        );
        let ess = ql_vol::evaluate_mingone_global_arbitrage_free_essvi(
            &ql_vol::MingoneGlobalArbitrageFreeEssviSpec {
                expiries: vec![1.0, 2.0],
                raw_theta_increments: vec![vol, vol],
                raw_rho_params: vec![-0.2, -0.1],
                raw_psi_fractions: vec![0.2, 0.3],
                test_log_moneyness_grid: vec![-1.0, 0.0, 1.0],
            },
        );
        let cb = ql_fx::price_cross_currency_convertible_and_warrant_bond(
            &ql_fx::CrossCurrencyConvertibleWarrantSpec {
                domestic_stock_spot: spot,
                domestic_conversion_strike: spot,
                spot_fx_for_per_dom: 1.1,
                fixed_contractual_fx_for_per_dom: 1.1,
                foreign_bond_par: 100.0,
                foreign_coupon_rate: 0.03,
                r_domestic: r,
                r_foreign: r - 0.005,
                stock_dividend_yield: 0.01,
                credit_spread: credit_bump + 0.02,
                stock_vol: vol,
                fx_vol: vol * 0.7,
                corr_stock_fx: 0.2,
                maturity: 5.0,
            },
        );
        let bank = ql_fixed_income::analyze_acharya_feds_bank_run_and_optimal_coco(
            &ql_fixed_income::AcharyaFedsBankRunAndCocoSpec {
                initial_securities_assets: 100.0,
                asset_effective_duration: 5.0,
                insured_deposits: 40.0,
                uninsured_deposits: 45.0,
                deposit_rate_beta: 0.5,
                deposit_operating_cost_rate: 0.005,
                at1_coco_principal: 5.0,
                initial_short_rate: r,
                parallel_rate_shock: inf,
                run_sensitivity_steepness: 30.0,
                fire_sale_liquidation_haircut: 0.1,
                coco_write_down_efficiency: 0.9,
            },
        );
        let xva = ql_xva::compute_albanese_crepey_mpor_and_balance_sheet_xva(
            &ql_xva::AlbaneseCrepeyMporBalanceSheetSpec {
                exposure_times: vec![1.0, 2.0],
                clean_mtm_profile: vec![1.0, 2.0],
                portfolio_vol_profile: vec![vol, vol],
                scheduled_cashflows_in_mpor: vec![0.01, 0.02],
                mpor_years: 0.05,
                csa_threshold: 0.1,
                risk_free_rate: r,
                cpty_hazard_rate: 0.02,
                cpty_lgd: 0.6,
                bank_hazard_rate: 0.01,
                bank_lgd: 0.6,
                economic_capital_profile: vec![0.05, 0.05],
                hurdle_rate: 0.1,
                initial_margin_profile: vec![0.1, 0.1],
                im_funding_spread: 0.01,
            },
        );
        let dh = ql_models::evaluate_double_heston_cumulants_and_kernel_slv(
            &ql_models::DoubleHestonKernelSlvSpec {
                spot,
                strike: spot,
                maturity: 1.0,
                risk_free_rate: r,
                dividend_yield: 0.0,
                v1_0: vol * vol * 0.6,
                kappa1: 2.0,
                theta1: vol * vol * 0.6,
                sigma1: 0.3,
                rho1: -0.6,
                v2_0: vol * vol * 0.4,
                kappa2: 0.5,
                theta2: vol * vol * 0.4,
                sigma2: 0.2,
                rho2: -0.3,
                target_dupire_local_vol: vol,
            },
        );
        let asian = ql_exotics::price_curran_rogers_shi_arithmetic_asian(
            &ql_exotics::CurranRogersShiArithmeticAsianSpec {
                spot,
                strike: spot,
                risk_free_rate: r,
                dividend_yield: 0.0,
                vol,
                fixing_times: vec![0.25, 0.5, 0.75, 1.0],
            },
        );
        let storage = ql_commodities::evaluate_warin_baeuerle_storage_rolling_hedge(
            &ql_commodities::WarinBaeuerleStorageRollingHedgeSpec {
                forward_curve: vec![spot, spot * 1.01, spot * 1.02],
                decision_times: vec![0.5, 1.0, 1.5],
                max_capacity: 10.0,
                initial_inventory: 2.0,
                max_injection_step: 1.0,
                max_withdrawal_step: 1.0,
                injection_loss_fraction: 0.02,
                spot_ou_volatility: vol,
                spot_ou_mean_reversion: 1.0,
                discount_rate: r,
            },
        );
        let micro = ql_portfolio::analyze_harris_roll_corwin_schultz_microstructure(
            &ql_portfolio::HarrisRollCorwinSchultzMicrostructureSpec {
                trade_prices: vec![spot, spot * 1.001, spot * 0.999, spot * 1.002],
                daily_highs: vec![spot * 1.01, spot * 1.02, spot * 1.015],
                daily_lows: vec![spot * 0.99, spot * 0.98, spot * 0.985],
                daily_dollar_volumes: vec![1e6, 1.1e6, 0.9e6],
                decision_price: spot,
                arrival_price: spot,
                child_fill_quantities: vec![10.0, 5.0],
                child_fill_prices: vec![spot * 1.001, spot * 1.002],
                total_order_quantity: 20.0,
                post_trade_benchmark_price: spot * 1.003,
                explicit_commissions: 0.1,
            },
        );
        let kelly = ql_portfolio::solve_busseti_boyd_risk_constrained_kelly(
            &ql_portfolio::BussetiBoydRiskConstrainedKellySpec {
                excess_return_means: vec![0.05 + spot * 1e-4 + inf * 0.1, 0.03 + r * 0.1],
                covariance: vec![
                    vec![vol * vol, vol * vol * 0.25],
                    vec![vol * vol * 0.25, (vol * vol * 0.75).max(1e-6)],
                ],
                risk_free_rate: r,
                drawdown_floor_alpha: 0.8,
                max_drawdown_prob_beta: 0.5,
                specialist_immediate_capital: 100.0,
                slow_moving_external_capital: 100.0,
                capital_arrival_rate: 0.1,
                dislocation_horizon_years: 1.0,
                forced_liquidation_notional: 10.0,
            },
        );
        let qp = ql_dml::solve_agrawal_boyd_differentiable_qp_layer(
            &ql_dml::AgrawalBoydDifferentiableQpSpec {
                q_matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
                q_linear: vec![-spot * 0.001, -vol],
                a_eq: vec![vec![1.0, 1.0]],
                b_eq: vec![1.0],
                g_ineq: vec![vec![-1.0, 0.0]],
                h_ineq: vec![0.0],
                barrier_tau: 1e-4,
            },
        );
        ftpl.equilibrium_price_level_p0.val()
            + g2.payer_swaption_pv.val()
            + ess.calibrated_thetas.iter().map(|x| x.val()).sum::<f64>()
            + cb.embedded_cross_currency_cb_price.val() * 1e-2
            + bank.economic_equity_post_coco_conversion.val() * 1e-2
            + xva.sticky_mpor_cva.val()
            + dh.double_heston_atm_vol.val()
            + asian.geometric_asian_call_price.val()
            + storage.stochastic_total_storage_value.val() * 1e-2
            + micro.perold_total_implementation_shortfall.val() * 1e-2
            + kelly.constrained_portfolio_volatility.val()
            + qp.optimal_objective_value.val()
            + inf
    }

    pub fn romer_taylor_damodaran_crp_christoffersen_var_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        Self::romer_taylor_damodaran_crp_christoffersen_var_live_metric_from_market_with_bump(
            m,
            credit_bump,
        )
    }
    pub fn romer_taylor_damodaran_crp_christoffersen_var_live_metric_from_market_with_bump(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03);
        let spot = m.eq_spot.max(1e-6);
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, x)| *x).unwrap_or(0.023);
        let taylor = ql_curves::evaluate_romer_veronesi_taylor_rule_macro_yield_curve(
            &ql_curves::RomerVeronesiTaylorRuleMacroCurveSpec {
                natural_real_rate_r_star: r - inf,
                inflation_target_pi_star: 0.02,
                current_inflation: inf,
                current_output_gap: 0.01,
                taylor_inflation_weight_phi_pi: 1.5,
                taylor_output_weight_phi_y: 0.5,
                inflation_mean_reversion: 0.4,
                output_gap_mean_reversion: 0.3,
                inflation_vol: vol,
                output_gap_vol: vol * 0.5,
                market_price_of_inflation_risk: 0.2,
                zlb_floor_rate: 0.0,
                maturities: vec![1.0, 5.0, 10.0],
            },
        );
        let dam = ql_credit::evaluate_damodaran_rsm_country_risk_and_debt_modification(
            &ql_credit::DamodaranRsmCountryRiskAndDebtModSpec {
                mature_market_erp: 0.05,
                sovereign_cds_spread: 0.012 + credit_bump,
                sovereign_default_drag_bps: 20.0,
                emerging_equity_vol: vol * 1.4,
                emerging_govt_bond_vol: vol.max(0.01),
                firm_domestic_revenue_share: 0.6,
                ebit_operating_income: 12.0,
                annual_interest_expense: 3.0,
                original_debt_carrying_value: 100.0,
                original_effective_interest_rate: r,
                modified_cashflows: vec![4.0, 4.0, 104.0],
                modified_payment_times: vec![1.0, 2.0, 3.0],
                modification_fees_paid: 1.0,
            },
        );
        let var = ql_portfolio::evaluate_christoffersen_bouchaud_var_backtest_and_cumulant_es(
            &ql_portfolio::ChristoffersenBouchaudVarBacktestAndCumulantSpec {
                total_days_t: 250,
                count_n00: 220.0,
                count_n01: 8.0,
                count_n10: 8.0,
                count_n11: 14.0,
                var_tail_probability_p: 0.05,
                portfolio_mean_loss: 0.0,
                portfolio_std_loss: 1.0 + vol,
                loss_skewness_gamma1: 0.2,
                loss_excess_kurtosis_gamma2: 0.4,
            },
        );
        let lvar = ql_portfolio::evaluate_foucault_vandermerwe_lob_schedule_and_lvar(
            &ql_portfolio::FoucaultVanDerMerweLobLiquidityLvarSpec {
                mid_price: spot,
                fundamental_volatility: vol,
                informed_trader_share_pi: 0.25,
                exponential_order_size_rate_lambda: 0.8,
                order_processing_fee_c: 0.0002,
                position_notional: 1_000_000.0,
                mean_relative_bid_ask_spread: 0.001,
                std_relative_bid_ask_spread: 0.0004,
                spread_tail_scaling_a: 1.5,
                confidence_z_alpha: 2.33,
                liquidation_horizon_days: 5.0,
                market_adv_notional: 10_000_000.0,
                permanent_impact_gamma: 0.02,
            },
        );
        let pt = ql_portfolio::evaluate_barberis_jin_wang_prospect_overhang_and_skew_anomaly(
            &ql_portfolio::BarberisJinWangProspectAnomalySpec {
                current_price: spot,
                reference_purchase_price: 95.0,
                return_volatility: vol,
                return_skewness: 0.8,
                prospect_loss_aversion_lambda: 2.25,
                prospect_diminishing_sensitivity_alpha: 0.88,
                prelec_probability_weighting_gamma: 0.7,
                mental_accounting_integration_b0: 0.4,
                capm_beta: 1.1,
                equity_risk_premium: 0.05,
            },
        );
        let regret = ql_portfolio::evaluate_hart_mascolell_regret_matching_allocation(
            &ql_portfolio::HartMasColellRegretMatchingPortfolioSpec {
                strategy_payoff_history: vec![
                    [0.02, 0.01, 0.03],
                    [0.01, 0.02, 0.0],
                    [0.03, 0.01, 0.02],
                ],
                realized_action_weights_history: vec![
                    [0.4, 0.3, 0.3],
                    [0.3, 0.4, 0.3],
                    [0.3, 0.3, 0.4],
                ],
                regret_inertia_mu: 1.0,
                exploration_epsilon: 0.01,
            },
        );
        let forensic = ql_ad::evaluate_beneish_wahlen_forensic_mscore_and_hkas32_bifurcation(
            &ql_ad::BeneishWahlenHkas32ForensicBifurcationSpec {
                dsri: 1.1,
                gmi: 1.05,
                aqi: 1.2,
                sgi: 1.15,
                depi: 1.03,
                sgai: 1.02,
                lvgi: 1.1,
                tata: 0.03,
                net_income: 10.0,
                operating_cash_flow: 8.0,
                average_total_assets: 200.0,
                convertible_issue_proceeds: 100.0,
                convertible_face_value: 100.0,
                contractual_coupon_rate: 0.03,
                comparable_straight_debt_yield: r + credit_bump,
                direct_issuance_transaction_costs: 1.0,
                maturity_years: 5,
            },
        );
        let camp = ql_portfolio::evaluate_colin_bacon_campisi_fixed_income_attribution(
            &ql_portfolio::ColinBaconCampisiFixedIncomeAttributionSpec {
                portfolio_total_return: 0.042 + credit_bump,
                benchmark_total_return: 0.035,
                portfolio_coupon_yield: 0.04,
                benchmark_coupon_yield: 0.035,
                portfolio_eff_duration: 5.2,
                benchmark_eff_duration: 5.0,
                portfolio_spread_duration: 4.0,
                benchmark_spread_duration: 4.2,
                parallel_treasury_shift: 0.005,
                portfolio_key_rate_twist_shift: 0.001,
                benchmark_key_rate_twist_shift: 0.001,
                portfolio_credit_spread_change: 0.002,
                benchmark_credit_spread_change: 0.002,
                fx_spot_return: 0.01,
                domestic_cash_rate: r,
                foreign_cash_rate: 0.032,
            },
        );
        let heck = ql_dml::evaluate_greene_wooldridge_heckman_survivorship_and_gmm(
            &ql_dml::GreeneWooldridgeHeckmanSelectionGmmSpec {
                selection_z_scores: vec![-0.5, 0.0, 0.5, 1.0],
                factor_exposures_x: vec![0.5, 1.0, 1.5, 2.0],
                observed_returns_y: vec![0.02, 0.04, 0.06, 0.09],
                gmm_moment_conditions_bar_g: [0.01, -0.02, 0.015],
                gmm_weighting_matrix_diag: [2.0, 1.5, 1.0],
                sample_size_n: 100.0,
            },
        );
        let en = ql_dml::solve_hastie_tibshirani_elastic_net_and_gcv(
            &ql_dml::HastieTibshiraniElasticNetGcvSpec {
                design_matrix_x: vec![
                    [1.0, 0.2, 0.1],
                    [1.0, 0.4, 0.5],
                    [1.0, 0.6, 0.2],
                    [1.0, 0.8, 0.9],
                ],
                response_y: vec![0.1, 0.2, 0.3, 0.4],
                l1_penalty_lambda1: 0.05 + credit_bump,
                l2_penalty_lambda2: 0.1,
                smoothing_eps: 1e-4,
            },
        );
        let cb = ql_fixed_income::evaluate_an_daoquan_clifford_chance_resettable_convertible(
            &ql_fixed_income::AnDaoquanCliffordChanceResettableCbSpec {
                spot_price: spot,
                initial_conversion_price: 100.0,
                book_value_per_share_bps: 80.0,
                avg_20day_price: spot * 0.95,
                face_value: 100.0,
                coupon_rate: 0.02,
                redemption_maturity_price: 105.0,
                put_back_trigger_ratio: 0.7,
                put_back_price: 103.0,
                soft_call_trigger_ratio: 1.3,
                soft_call_price: 105.0,
                stock_volatility: vol,
                risk_free_rate: r,
                credit_spread: 0.01 + credit_bump,
                maturity: 3.0,
            },
        );
        let bail = ql_credit::evaluate_fsb_finma_gsib_resolution_bail_in_waterfall(
            &ql_credit::FsbFinmaGsibResolutionBailInWaterfallSpec {
                total_rwa: 1000.0,
                pre_resolution_cet1: 120.0,
                at1_principal: 50.0,
                tier2_principal: 60.0,
                senior_tlac_principal: 100.0,
                valuation_loss_shock: 0.08 + credit_bump,
                target_post_resolution_cet1_ratio: 0.1,
                liquidation_insolvency_haircut: 0.4,
                existing_share_count: 100.0,
                bail_in_conversion_discount: 0.2 + credit_bump * 0.01,
            },
        );
        taylor.macro_affine_zero_yields.iter().copied().sum::<f64>()
            + dam.country_risk_premium_crp
            + var.cornish_fisher_var
            + lvar.total_liquidity_adjusted_var_lvar * 1e-6
            + pt.bjw_equilibrium_expected_excess_return
            + regret.hannan_external_regret_per_period
            + forensic.beneish_8var_m_score
            + camp.total_active_local_return
            + heck.heckman_corrected_alpha
            + en.generalized_cross_validation_gcv
            + cb.three_line_duplex_fair_value * 1e-3
            + bail.post_bail_in_creditor_ownership_share
            + bail.ncwo_compensation_shortfall * 1e-3
            + dam.total_cost_of_equity_ke
            + dam.extinguishment_gain_or_loss * 1e-3
            + forensic.manipulation_probit_probability
            + forensic.year1_effective_interest_expense * 1e-2
            + camp.active_security_selection_residual
            + camp.karnosky_singer_unhedged_currency_surprise
    }
    pub fn romer_taylor_damodaran_crp_christoffersen_var_live_metric(&self) -> f64 {
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::romer_taylor_damodaran_crp_christoffersen_var_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03);
        let inf = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let f = m.eq_spot.max(1.0);
        let vol = m.eq_vol.max(0.1);
        let c = ql_models::evaluate_gairat_polala_cheyette_local_vol_caplet(
            &ql_models::GairatPolalaCheyetteCapletSpec {
                forward_rate: r,
                strike_rate: r + credit_bump,
                fixing_time_t1: 1.0,
                payment_time_t2: 1.5,
                mean_reversion_kappa: 0.1,
                initial_sigma0: vol * 0.1,
                initial_skew_b: -0.2,
                target_atm_normal_vol: vol * 0.01,
                target_otm_normal_vol: vol * 0.011,
                discount_factor_t2: (-r * 1.5).exp(),
            },
        );
        let q = ql_models::evaluate_gatheral_jusselin_quadratic_rough_heston(
            &ql_models::GatheralJusselinQuadraticRoughHestonSpec {
                forward: f,
                strike: f * 1.01,
                expiry: 1.0,
                hurst_h: 0.1,
                weight_a: 0.2,
                asymmetry_center_b: 0.3,
                floor_variance_c: 0.02,
                initial_state_z0: 0.1,
                kernel_eta: 0.5,
                mean_reversion_lambda: 1.0,
            },
        );
        let z = ql_vol::evaluate_andreasen_huge_one_step_zabr_calibration(
            &ql_vol::AndreasenHugeOneStepZabrSpec {
                forward: f,
                strike: f * 1.01,
                expiry: 1.0,
                alpha: vol * 0.2,
                beta: 0.5,
                rho: -0.3,
                nu: 0.4,
                discount_factor: 1.0,
            },
        );
        let n = ql_curves::evaluate_adb_rbi_ndirs_and_asset_swap_basis(
            &ql_curves::AdbRbiNdirsAssetSwapSpec {
                payment_tenors: vec![1.0, 2.0],
                onshore_forward_rates: vec![r, r + 0.002],
                usd_ois_discount_factors: vec![1.0, 0.95],
                local_rate_vol: 0.1,
                ndf_fx_vol: vol * 0.1,
                rate_fx_corr: 0.2,
                sovereign_bond_yield: r + inf,
                usd_sofr_yield: r,
                ndf_ccs_basis_spread: 0.005,
                withholding_tax_rate: 0.1,
            },
        );
        let x = ql_xva::evaluate_turlakov_green_kenyon_wwr_and_kva(
            &ql_xva::TurlakovGreenKenyonWwrKvaSpec {
                time_grid: vec![1.0, 2.0],
                independent_epe_profile: vec![f * 0.01, f * 0.012],
                pfe_addon_profile: vec![f * 0.005, f * 0.006],
                base_cpty_hazard: 0.02 + credit_bump,
                hazard_log_vol: 0.3,
                exposure_vol: 0.2,
                wwr_correlation_rho: 0.5,
                lgd: 0.6,
                risk_free_rate: r,
                bank_hazard_rate: 0.01,
                ccr_risk_weight: 0.08,
                shareholder_cost_of_capital: 0.1,
                capital_return_offset_beta: 0.5,
            },
        );
        let e = ql_portfolio::evaluate_ahmadi_javid_entropic_var_portfolio(
            &ql_portfolio::AhmadiJavidEntropicVarPortfolioSpec {
                scenario_returns: vec![vec![0.01, -0.02], vec![-0.03, 0.01]],
                scenario_probabilities: vec![0.5, 0.5],
                initial_weights: vec![0.6, 0.4],
                tail_alpha: 0.05,
                risk_aversion_lambda: 1.0,
            },
        );
        let h = ql_portfolio::evaluate_cotton_schur_hrp_and_rmt_nco(
            &ql_portfolio::CottonSchurHrpRmtNcoSpec {
                covariance_4x4: [
                    [0.04, 0.01, 0.0, 0.0],
                    [0.01, 0.05, 0.0, 0.0],
                    [0.0, 0.0, 0.03, 0.01],
                    [0.0, 0.0, 0.01, 0.06],
                ],
                schur_gamma: 0.5,
                sample_aspect_ratio_q: 0.5,
            },
        );
        let b = ql_portfolio::evaluate_tsay_bekk_garch_volatility_spillover(
            &ql_portfolio::TsayBekkGarchSpilloverSpec {
                chol_c: [[0.1, 0.0], [0.02, 0.1]],
                arch_a: [[0.1, 0.02], [0.01, 0.1]],
                garch_b: [[0.8, 0.0], [0.0, 0.7]],
                initial_cov: [[0.04, 0.0], [0.0, 0.04]],
                return_innovations: vec![[0.01, -0.02], [0.02, 0.01]],
            },
        );
        let ql = ql_fixed_income::evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(
            &ql_fixed_income::MayerBrownSsgaCappedCallMakeWholeSpec {
                bond_par: 100.0,
                spot_price: f,
                base_conversion_price: f * 0.9,
                cap_strike_price: f * 1.2,
                volatility: vol,
                risk_free_rate: r,
                maturity_years: 5.0,
                acquisition_offer_price: f * 1.1,
                remaining_years_at_acquisition: 2.0,
                tax_rate_section_163j: 0.21,
            },
        );
        let hhw = ql_pde::evaluate_haentjens_inthout_heston_hull_white_adi(
            &ql_pde::HaentjensInthoutHhwAdiSpec {
                spot: f,
                strike: f,
                expiry: 1.0,
                initial_var_v0: 0.04,
                mean_rev_kappa_v: 1.0,
                long_var_theta_v: 0.04,
                vol_of_vol_sigma_v: 0.3,
                initial_rate_r0: r,
                rate_mean_rev_kappa_r: 0.1,
                rate_vol_sigma_r: 0.01,
                corr_spot_vol: 0.2,
                corr_spot_rate: 0.3,
            },
        );
        let qual = ql_models::evaluate_schilit_sloan_dechow_dichev_earnings_quality(
            &ql_models::SchilitSloanDechowDichevQualitySpec {
                revenue_0: 100.0,
                revenue_1: 105.0,
                cogs_0: 60.0,
                cogs_1: 63.0,
                accounts_receivable_0: 10.0,
                accounts_receivable_1: 12.0,
                inventory_0: 8.0,
                inventory_1: 9.0,
                accounts_payable_0: 7.0,
                accounts_payable_1: 8.0,
                net_operating_assets_0: 50.0,
                net_operating_assets_1: 52.0,
                net_income_1: 8.0,
                cfo_prev: 7.0,
                cfo_curr: 7.5,
                cfo_next: 8.0,
            },
        );
        let cr = ql_credit::evaluate_brigo_morini_cmcds_and_credit_index_option(
            &ql_credit::BrigoMoriniCmcdsCreditOptionSpec {
                reset_times: vec![1.0, 2.0],
                forward_cds_spreads: vec![0.01, 0.012],
                spread_volatility: 0.3,
                annuity_convexity_gamma: 0.5,
                flat_par_cds_spread: 0.011,
                option_expiry: 1.0,
                option_strike_spread: 0.01,
                risky_pv01_annuity: 4.0,
                hazard_rate: 0.02,
                recovery_rate: 0.4,
                risk_free_rate: r,
            },
        );
        c.deterministic_auxiliary_y_bar.val()
            + q.qrh_call_price.val()
            + z.one_step_fd_call_price.val()
            + n.usd_asset_swap_pickup_bps.val() * 1e-4
            + x.total_wwr_cva_plus_kva.val()
            + e.exact_chernoff_evar.val()
            + h.schur_hrp_portfolio_variance.val()
            + b.terminal_var_asset1.val()
            + ql.capped_call_overlay_gross_cost.val() * 1e-2
            + cr.morini_brigo_payer_swaption_pv.val()
            + hhw.hhw_adi_call_price.val()
            + qual.schilit_shenanigan_composite_z.val()
    }
    pub fn tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.10);
        let repo = ql_curves::evaluate_tuckman_repo_specialness_and_oftr_basis(
            &ql_curves::TuckmanRepoSpecialnessBasisSpec {
                otr_dirty_price: spot,
                otr_yield: r,
                oftr_yield: r + 0.0015,
                modified_duration: 7.5,
                gc_repo_rate: r + 0.01,
                special_repo_rate: r - 0.005,
                specialness_horizon_years: 0.25,
                repo_haircut_fraction: 0.02,
                equity_funding_rate: r + 0.02,
            },
        );
        let bk = ql_credit::evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(
            &ql_credit::CrepeyBurgardKjaerBilateralFvaSpec {
                time_grid: vec![0.5, 1.0, 2.0],
                epe_profile: vec![10.0, 12.0, 8.0],
                ene_profile: vec![6.0, 7.0, 5.0],
                cpty_hazard_rate: 0.02 + credit_bump,
                bank_hazard_rate: 0.015 + inflation * 0.1,
                cpty_recovery_rate: 0.4,
                bank_recovery_rate: 0.4,
                risk_free_rate: r,
                liquidity_funding_spread: 0.008,
                burgard_kjaer_alpha: 0.5,
            },
        );
        let decs = ql_fixed_income::evaluate_despiegeleer_mandatory_convertible_decs_and_percs(
            &ql_fixed_income::DeSpiegeleerMandatoryConvertibleDecsSpec {
                spot_price: spot,
                par_amount: 100.0,
                lower_strike_kl: 90.0,
                upper_strike_ku: 120.0,
                decs_coupon_rate: 0.06,
                risk_free_rate: r,
                dividend_yield: 0.02,
                credit_spread: 0.015 + credit_bump,
                lower_strike_vol: vol + 0.03,
                upper_strike_vol: vol,
                maturity_years: 3.0,
            },
        );
        let mbs = ql_fixed_income::evaluate_davidson_levin_four_factor_prepayment_io_po(
            &ql_fixed_income::DavidsonLevinFourFactorIoPoSpec {
                initial_balance: 100.0,
                gross_wac_rate: 0.06,
                net_passthrough_coupon: 0.055,
                current_mortgage_rate: r + 0.02,
                base_turnover_cpr: 0.06,
                max_refi_cpr: 0.35,
                burnout_sensitivity_beta: 0.5,
                seasonality_amplitude: 0.1,
                oas_discount_rate: r + 0.01,
                num_periods: 12,
            },
        );
        let kyle = ql_portfolio::evaluate_foucault_kyle_lambda_and_ekop_pin(
            &ql_portfolio::FoucaultKyleEkopPinSpec {
                prior_price_p0: spot,
                insider_valuation_v: spot * 1.05,
                fundamental_std_sigma_v: 2.0 + vol,
                noise_order_std_sigma_u: 10.0,
                event_prob_alpha: 0.3,
                bad_news_prob_delta: 0.4,
                informed_arrival_mu: 40.0,
                uninformed_buy_eps_b: 50.0,
                uninformed_sell_eps_s: 50.0,
                observed_buys_b: 75.0 + r * 10.0,
                observed_sells_s: 48.0,
            },
        );
        let gp = ql_portfolio::evaluate_garleanu_pedersen_dynamic_aim_portfolio(
            &ql_portfolio::GarleanuPedersenDynamicAimSpec {
                current_weights: [0.2, 0.1],
                fast_alpha_signal: [0.04 + r * 0.1, 0.02],
                slow_alpha_signal: [0.03, 0.025],
                fast_decay_phi: 0.5,
                slow_decay_phi: 0.1,
                covariance_2x2: [[0.04, 0.01], [0.01, 0.05]],
                risk_aversion_gamma: 2.0,
                transaction_cost_lambda: 0.5,
                discount_rate_rho: 0.05,
            },
        );
        let gjr = ql_portfolio::evaluate_tsay_gjr_and_dcc_garch_dynamic_correlation(
            &ql_portfolio::TsayGjrDccGarchSpec {
                bivariate_returns: vec![
                    [-0.02, -0.015],
                    [0.012 + r * 0.05, -0.01],
                    [-0.018, -0.022],
                ],
                omega: [0.00002, 0.00003],
                alpha: [0.05, 0.06],
                gjr_gamma: [0.08, 0.07],
                beta: [0.85, 0.84],
                unconditional_corr: 0.5,
                dcc_a: 0.05,
                dcc_b: 0.90,
            },
        );
        let lppl = ql_portfolio::evaluate_sornette_jls_lppl_crash_hazard(
            &ql_portfolio::SornetteJlsLpplCrashHazardSpec {
                current_time_t: 1.0,
                critical_time_tc: 1.5,
                horizon_dt: 0.1,
                peak_log_price_a: 5.0 + r * 0.1,
                power_law_slope_b: -0.4,
                log_periodic_amp_c: 0.03,
                critical_exponent_m: 0.5,
                log_periodic_omega: 6.5,
                phase_phi: 1.0,
                crash_loss_kappa: 0.2,
            },
        );
        let spec = ql_portfolio::evaluate_follmer_schied_entropic_and_spectral_risk(
            &ql_portfolio::FollmerSchiedEntropicSpectralRiskSpec {
                ordered_losses: vec![0.5, 1.2, 2.5, 4.0 + r],
                entropic_theta: 0.5,
                spectral_power_gamma: 2.5,
            },
        );
        let brin = ql_portfolio::evaluate_bacon_brinson_fachler_and_carino_linking(
            &ql_portfolio::BaconBrinsonFachlerCarinoSpec {
                portfolio_sector_weights: vec![0.6, 0.4],
                benchmark_sector_weights: vec![0.5, 0.5],
                portfolio_sector_returns: vec![0.08 + r * 0.1, 0.03],
                benchmark_sector_returns: vec![0.06, 0.02],
                multi_period_port_returns: vec![0.04 + r * 0.1, 0.03, -0.01],
                multi_period_bench_returns: vec![0.03, 0.02, -0.015],
            },
        );
        let prel = ql_portfolio::evaluate_wakker_prelec_rdu_and_choquet_ambiguity(
            &ql_portfolio::WakkerPrelecChoquetProspectSpec {
                sorted_outcomes: [-10.0, -2.0, 8.0, 20.0 + r * 10.0],
                state_probabilities: [0.1, 0.3, 0.4, 0.2],
                reference_point: 0.0,
                prelec_curvature_alpha: 0.65,
                prelec_elevation_beta: 1.0,
                tk_value_exponent: 0.88,
                loss_aversion_lambda: 2.25,
                ambiguity_aversion_delta: 0.15,
            },
        );
        let sig = ql_dml::evaluate_bayer_path_signature_chen_and_levy_area(
            &ql_dml::BayerPathSignatureChenLevyAreaSpec {
                time_grid: vec![0.0, 0.5, 1.0],
                spot_path: vec![spot, spot * 1.04, spot * (1.02 + r * 0.1)],
                linear_functional_weights: [0.1, 0.2, 0.3, 0.4, -0.2, 0.1],
            },
        );
        repo.structural_liquidity_yield_spread_bps * 0.01
            + bk.total_burgard_kjaer_xva_charge
            + decs.decs_mandatory_fair_value * 0.01
            + mbs.passthrough_present_value * 0.01
            + kyle.ekop_pin_probability
            + gp.expected_net_utility_gain
            + gjr.terminal_equal_weight_portfolio_vol
            + lppl.lppl_fitted_price * 0.01
            + spec.entropic_convex_risk_measure * 0.1
            + brin.carino_smoothed_active_return
            + prel.prospect_certainty_equivalent * 0.1
            + sig.signature_linear_payoff * 0.01
    }
    pub fn tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);
        let hw = ql_curves::evaluate_hagan_west_monotone_convex_forward_curve(
            &ql_curves::HaganWestMonotoneConvexCurveSpec {
                knot_tenors: [0.0, 1.0, 2.0, 3.0],
                discrete_forwards: [r, r + 0.008, r + 0.012],
                eval_tenor: 1.4,
            },
        );
        let ctd = ql_fixed_income::evaluate_fabozzi_treasury_futures_ctd_and_switch_option(
            &ql_fixed_income::FabozziTreasuryFuturesCtdSwitchSpec {
                bond_coupons: [0.06, 0.045],
                bond_maturities: [7.0, 9.5],
                bond_clean_prices: [101.0 - r * 20.0, 94.5 - r * 25.0],
                bond_modified_durations: [5.8, 7.6],
                futures_price: 100.2 - r * 15.0,
                delivery_horizon_years: 0.25,
                repo_rate: r,
                yield_volatility: 0.01,
            },
        );
        let afv = ql_fixed_income::evaluate_ayache_forsyth_vetzal_default_coupled_cb(
            &ql_fixed_income::AyacheForsythVetzalConvertibleSpec {
                spot_price: spot,
                bond_par: 100.0,
                coupon_rate: 0.04,
                conversion_ratio: 1.0,
                maturity_years: 3.0,
                risk_free_rate: r,
                dividend_yield: 0.015,
                equity_volatility: vol,
                hazard_rate: 0.025 + credit_bump,
                bond_recovery_rate: 0.4,
                equity_default_jump_loss_eta: 0.7,
            },
        );
        let kth = ql_credit::evaluate_lando_schonbucher_kth_to_default_basket(
            &ql_credit::LandoSchonbucherBasketDefaultSpec {
                obligor_hazard_rates: [
                    0.015 + credit_bump,
                    0.02 + credit_bump,
                    0.025 + credit_bump,
                    0.03 + credit_bump,
                ],
                factor_loading_beta: 0.5,
                recovery_rate: 0.4,
                risk_free_rate: r,
                maturity_years: 5.0,
            },
        );
        let sacva = ql_xva::evaluate_gregory_sa_cva_and_euler_netting_allocation(
            &ql_xva::GregorySaCvaEulerNettingSpec {
                cva_delta_sensitivities: [10.0 + r * 10.0, 8.0],
                cva_hedge_sensitivities: [6.0, 5.0],
                delta_risk_weights: [0.05, 0.06],
                intra_bucket_corr_rho: 0.5,
                cva_vega_sensitivity: 4.0 + vol * 5.0,
                vega_risk_weight: 0.1,
                subportfolio_exposures: [20.0, 15.0, 10.0],
                exposure_corr_3x3: [[1.0, 0.3, 0.1], [0.3, 1.0, 0.2], [0.1, 0.2, 1.0]],
                cpty_hazard_rate: 0.02 + credit_bump,
                cpty_lgd: 0.6,
                horizon_years: 3.0,
            },
        );
        let kirk = ql_commodities::evaluate_eydeland_gibson_schwartz_and_kirk_spread(
            &ql_commodities::EydelandGibsonSchwartzKirkSpreadSpec {
                spot_price: 55.0 + inflation * 50.0,
                initial_convenience_yield: 0.06,
                mean_reversion_kappa: 1.2,
                long_run_convenience_yield: 0.05,
                spot_volatility: 0.32,
                convenience_yield_vol: 0.20,
                spot_cy_correlation: 0.4,
                risk_free_rate: r,
                maturity_years: 1.0,
                fuel_futures_price: 6.0,
                heat_rate_multiplier: 7.5,
                fixed_strike_cost: 5.0,
                fuel_volatility: 0.28,
                power_fuel_correlation: 0.6,
            },
        );
        let fx = ql_fx::evaluate_clark_wystup_fx_vanna_volga_smile(
            &ql_fx::ClarkWystupFxVannaVolgaSmileSpec {
                spot_fx: 1.10 + r * 0.1,
                domestic_rate: r,
                foreign_rate: 0.015,
                expiry_years: 1.0,
                atm_vol: 0.10 + vol * 0.1,
                risk_reversal_25d: 0.012,
                butterfly_25d: 0.004,
                target_strike: 1.14,
            },
        );
        let nap = ql_exotics::evaluate_bouzoubaa_napoleon_and_reverse_cliquet(
            &ql_exotics::BouzoubaaNapoleonReverseCliquetSpec {
                period_returns: vec![0.03, -0.04 - credit_bump, 0.02, -0.025],
                napoleon_fixed_coupon: 0.08,
                global_floor: 0.0,
                local_cap: 0.04,
                local_floor: -0.04,
                forward_volatility: vol,
                forward_skew: -0.15,
                risk_free_rate: r,
                maturity_years: 2.0,
            },
        );
        let tc = ql_portfolio::evaluate_grinold_qian_transfer_coefficient_and_flam(
            &ql_portfolio::GrinoldQianTransferCoefficientFlAmSpec {
                information_coefficient_ic: 0.06,
                breadth_br: 100.0,
                target_active_risk: 0.04 + r * 0.1,
                residual_volatilities: [0.15, 0.20, 0.18],
                alpha_forecasts: [0.03, 0.02, -0.025],
                constrained_active_weights: [0.12, 0.05, -0.08],
            },
        );
        let ep = ql_portfolio::evaluate_meucci_relative_entropy_view_pooling(
            &ql_portfolio::MeucciEntropyPoolingViewSpec {
                prior_probabilities: [0.25, 0.25, 0.25, 0.25],
                scenario_returns: [-0.04, 0.01, 0.03, 0.08 + r * 0.1],
                target_view_return: 0.035,
                view_confidence: 0.75,
            },
        );
        let krr = ql_dml::evaluate_hastie_kernel_ridge_and_exact_loocv(
            &ql_dml::HastieKernelRidgeLoocvSpec {
                feature_x: [-1.0, -0.3, 0.4, 1.1],
                target_y: [0.8 + r, 0.2, 0.5, 1.3],
                rbf_bandwidth_ell: 0.8,
                ridge_penalty_lambda: 0.15,
            },
        );
        let eva = ql_models::evaluate_koller_mckinsey_kvd_dcf_and_eva_valuation(
            &ql_models::KollerMckinseyKvdEvaValuationSpec {
                invested_capital_0: 500.0,
                noplat_1: 75.0 + spot * 0.05,
                ronic: 0.16,
                wacc: 0.09 + r * 0.2,
                perpetual_growth_g: 0.04,
                net_debt: 150.0,
                diluted_shares: 50.0,
            },
        );
        hw.eval_discount_factor
            + ctd.switch_adjusted_theoretical_futures * 0.01
            + afv.afv_convertible_fair_value * 0.01
            + kth.first_to_default_spread_bps * 0.001
            + sacva.total_sa_cva_capital
            + kirk.kirk_spark_spread_call_pv * 0.1
            + fx.vanna_volga_call_price * 10.0
            + nap.napoleon_capital_protected_note_pv
            + tc.expected_constrained_alpha * 10.0
            + ep.posterior_expected_return * 10.0
            + krr.loocv_shortcut_mse
            + eva.intrinsic_share_price * 0.1
    }
    pub fn hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let csa = ql_curves::evaluate_brigo_gregory_multicurrency_csa_ctd(
            &ql_curves::BrigoGregoryMultiCurrencyCsaCtdSpec {
                domestic_ois_rate: r,
                foreign_ois_rate: r - 0.012,
                xccy_basis_spread: -0.0025 - credit_bump * 0.2,
                csa_spread_volatility: 0.0065 + vol * 0.01,
                spot_fx: 1.085,
                maturity_years: 5.0,
                promised_cashflow: 100.0,
            },
        );
        let oas = ql_fixed_income::evaluate_fabozzi_callable_putable_oas_and_convexity(
            &ql_fixed_income::FabozziCallablePutableOasSpec {
                initial_short_rate: r,
                rate_drift_per_year: 0.002,
                short_rate_vol: 0.18 + vol * 0.1,
                oas_spread: 0.006 + credit_bump,
                coupon_rate: 0.055,
                call_strike_price: 100.5,
                put_strike_price: 99.5,
                maturity_years: 3.0,
                yield_shift_dy: 0.0025,
            },
        );
        let jlt = ql_credit::evaluate_lando_jlt_markov_rating_migration(
            &ql_credit::LandoJltMarkovMigrationSpec {
                ig_to_hy_rate: 0.045 + credit_bump * 2.0,
                ig_to_default_rate: 0.002 + credit_bump * 0.5,
                hy_to_ig_rate: 0.030,
                hy_to_ccc_rate: 0.080 + credit_bump * 2.0,
                hy_to_default_rate: 0.025 + credit_bump,
                ccc_to_hy_rate: 0.040,
                ccc_to_default_rate: 0.140 + credit_bump * 3.0,
                risk_premium_multiplier_pi: 1.45,
                recovery_rate: 0.40,
                risk_free_rate: r,
                horizon_years: 3.0,
            },
        );
        let corr_vs = ql_vol::evaluate_gatheral_bergomi_variance_swap_and_corridor(
            &ql_vol::GatheralBergomiVarianceSwapCorridorSpec {
                forward_price: spot,
                risk_free_rate: r,
                maturity_years: 1.0,
                atm_vol: vol,
                skew_slope_beta: -0.15,
                smile_curvature_gamma: 0.08,
                strike_grid: vec![spot * 0.80, spot * 0.90, spot, spot * 1.10, spot * 1.20],
                corridor_lower_barrier: spot * 0.85,
                corridor_upper_barrier: spot * 1.15,
            },
        );
        let lsm = ql_pde::evaluate_glasserman_lsm_and_rogers_dual_bound(
            &ql_pde::GlassermanLsmDualUpperBoundSpec {
                spot_price: spot,
                strike_price: 100.0,
                risk_free_rate: r,
                volatility: vol,
                maturity_years: 1.0,
                shock_paths_t1: [-1.35, -0.75, -0.25, 0.25, 0.75, 1.35],
                shock_paths_t2: [-0.95, 0.85, -0.45, 0.45, -0.85, 0.95],
            },
        );
        let lob = ql_portfolio::evaluate_cartea_jaimungal_lob_market_making(
            &ql_portfolio::CarteaJaimungalLobMarketMakingSpec {
                mid_price: spot,
                current_inventory_q: 3.0,
                running_inventory_penalty_phi: 0.02,
                terminal_inventory_penalty_alpha: 0.05,
                risk_aversion_gamma: 0.10,
                asset_volatility: vol,
                order_arrival_baseline_a: 12.0,
                lob_depth_decay_kappa: 1.5,
                remaining_horizon_tau: 0.5,
            },
        );
        let clm = ql_portfolio::evaluate_campbell_lo_mackinlay_vr_and_campbell_shiller(
            &ql_portfolio::CampbellLoMackinlayVrPresentValueSpec {
                returns: vec![
                    0.012,
                    0.015 - credit_bump,
                    -0.004,
                    0.009,
                    0.014,
                    -0.003,
                    0.011,
                    0.013,
                ],
                steady_state_dividend_yield: 0.032 + r * 0.1,
                expected_dividend_growth: 0.045 + inflation * 0.2,
                expected_discount_rate: 0.078 + r * 0.5,
                dividend_growth_ar1_phi: 0.35,
                discount_rate_ar1_phi: 0.85,
            },
        );
        let dsr = ql_dml::evaluate_lopez_de_prado_deflated_sharpe_and_purged_cv(
            &ql_dml::LopezDePradoDsrPurgedCvSpec {
                observed_sharpe_ratio: 1.45 - credit_bump * 5.0,
                sample_length_t: 60.0,
                return_skewness: -0.65,
                return_kurtosis: 4.5,
                num_independent_trials_k: 25.0,
                trials_sharpe_std: 0.45 + vol * 0.1,
                total_cv_samples_n: 1000.0,
                cv_folds_k: 5.0,
                label_overlap_horizon_h: 5.0,
                embargo_fraction_pct: 0.01,
            },
        );
        let erc = ql_portfolio::evaluate_roncalli_spinu_erc_risk_parity(
            &ql_portfolio::RoncalliSpinuErcRiskParitySpec {
                covariance_3x3: [
                    [vol * vol, 0.008, 0.002],
                    [0.008, 0.0225 + credit_bump * 0.1, 0.003],
                    [0.002, 0.003, 0.0064],
                ],
                target_risk_budgets: [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            },
        );
        let ilm = ql_portfolio::evaluate_ilmanen_carry_rolldown_and_vmc_factors(
            &ql_portfolio::IlmanenCarryRolldownVmcFactorSpec {
                current_yield_yt: r + 0.015,
                rolled_horizon_yield: r + 0.011,
                repo_financing_rate: r,
                modified_duration: 7.5,
                bond_convexity: 68.0,
                expected_yield_change: 0.001 + credit_bump * 0.2,
                yield_volatility: 0.0085 + vol * 0.01,
                value_style_sharpe: 0.45,
                momentum_style_sharpe: 0.50,
                carry_style_sharpe: 0.48,
                value_momentum_corr_rho: -0.45,
            },
        );
        let pen = ql_models::evaluate_penman_nissim_dupont_and_reoi_valuation(
            &ql_models::PenmanNissimDuPontReoiValuationSpec {
                net_operating_assets_0: 500.0,
                net_financial_obligations_0: 180.0,
                sales_revenue_1: 650.0 + spot * 0.5,
                after_tax_operating_income_1: 72.0 + inflation * 50.0,
                after_tax_net_financial_expense_1: 9.0 + credit_bump * 100.0,
                required_return_operations_rho: 0.09 + r * 0.2,
                perpetual_growth_g: 0.03,
            },
        );
        let mc = ql_models::evaluate_ljungqvist_sargent_mccall_bellman_search(
            &ql_models::LjungqvistSargentMccallSearchSpec {
                offer_states: [40.0, 55.0, 72.0, 95.0 + spot * 0.05],
                offer_probabilities: [0.25, 0.35, 0.25, 0.15],
                waiting_compensation_c: 32.0 + inflation * 20.0,
                discount_factor_beta: (1.0 / (1.0 + r)).clamp(0.80, 0.99),
            },
        );

        csa.ctd_csa_discounted_pv * 0.01
            + oas.callable_bond_fair_value * 0.01
            + jlt.ig_credit_spread_bps * 0.001
            + corr_vs.replicated_volatility_strike
            + lsm.lsm_control_variate_price * 0.01
            + lob.reservation_mid_price * 0.01
            + clm.equilibrium_log_dividend_price_ratio * 0.1
            + dsr.deflated_sharpe_ratio_dsr
            + erc.erc_portfolio_volatility
            + ilm.total_expected_excess_return
            + pen.intrinsic_common_equity_value * 0.001
            + mc.reservation_threshold_w_bar * 0.01
    }
    pub fn brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.05);
        let dx0 = (spot / 100.0 - 1.0) * 0.005 + credit_bump * 0.1;
        let dx1 = (r - 0.03) * 0.005 + credit_bump * 0.1;
        let dx2 = (vol - 0.20) * 0.005 - credit_bump * 0.1;
        let dx3 = credit_bump * 0.01 + inflation * 0.001;

        let nss = ql_curves::evaluate_brigo_tuckman_nss_and_key_rate_durations(
            &ql_curves::BrigoTuckmanNssKeyRateDurationSpec {
                beta0_level: 0.045 + dx0,
                beta1_slope: -0.015 + dx1,
                beta2_hump1: 0.018 + dx2,
                beta3_hump2: -0.008 + dx3,
                tau1_decay: 1.8,
                tau2_decay: 5.5,
                eval_maturity: 4.0,
                cashflow_weights_2y_5y_10y: [0.25, 0.35, 0.40],
            },
        );
        let tips = ql_fixed_income::evaluate_tuckman_fabozzi_tips_deflation_floor_and_bei(
            &ql_fixed_income::TuckmanFabozziTipsDeflationFloorSpec {
                base_cpi_index: 295.0,
                current_cpi_index: 300.0 + (spot - 100.0) * 0.1,
                nominal_treasury_yield: 0.042 + dx1,
                real_tips_yield: 0.019 + dx2,
                real_coupon_rate: 0.0175,
                cpi_volatility: 0.025 + dx3,
                inflation_risk_premium: 0.003,
                tips_liquidity_premium: 0.0015,
                maturity_years: 5.0,
            },
        );
        let rec = ql_credit::evaluate_schonbucher_lando_rmv_rt_rfv_recovery_triangle(
            &ql_credit::SchonbucherLandoRecoveryConventionSpec {
                par_amount: 100.0,
                risk_free_rate: 0.04 + dx0,
                hazard_rate_lambda: 0.035 + dx1,
                recovery_rate: 0.40 + dx2,
                coupon_rate: 0.05 + dx3,
                maturity_years: 5.0,
            },
        );
        let slv = ql_vol::evaluate_gatheral_bergomi_sticky_vs_local_vol_delta(
            &ql_vol::GatheralBergomiStickyLocalVolDeltaSpec {
                spot_price: spot,
                strike_price: 100.0,
                risk_free_rate: r,
                dividend_yield: 0.01,
                forward_start_t1: 0.5,
                maturity_t2: 1.5,
                atm_vol: vol,
                log_moneyness_skew_s0: -0.15,
                sv_mean_reversion_kappa: 1.2,
            },
        );
        let mva = ql_xva::evaluate_gregory_bilateral_simm_vs_ccp_mva(
            &ql_xva::GregoryBilateralSimmVsCcpMvaSpec {
                portfolio_annual_volatility: 15.0 + dx0 * 10.0,
                bilateral_mpor_years: 10.0 / 250.0,
                ccp_mpor_years: 5.0 / 250.0,
                exposure_decay_kappa: 0.12 + dx1,
                im_funding_spread: 0.0085 + dx2,
                cpty_hazard_rate: 0.015 + dx3,
                cpty_lgd: 0.60,
                risk_free_rate: r,
                ccp_default_fund_contribution: 0.35,
                capital_hurdle_rate: 0.10,
                maturity_years: 5.0,
            },
        );
        let ss = ql_commodities::evaluate_schwartz_smith_two_factor_and_kalman_update(
            &ql_commodities::SchwartzSmithTwoFactorKalmanSpec {
                prior_short_term_chi0: 0.08 + dx0,
                prior_equilibrium_xi0: 4.20 + dx1,
                mean_reversion_kappa: 1.4 + dx2,
                equilibrium_drift_mu: 0.01,
                short_term_vol_sigma_chi: 0.28 + dx3,
                equilibrium_vol_sigma_xi: 0.14,
                factor_correlation_rho: 0.30,
                maturity_years: 1.0,
                prior_var_p11: 0.04,
                prior_var_p22: 0.02,
                observed_log_futures: 4.32,
                measurement_noise_var_r: 0.005,
            },
        );
        let dnt = ql_fx::evaluate_wystup_clark_fx_double_no_touch_vanna_volga(
            &ql_fx::WystupClarkFxDoubleNoTouchSpec {
                spot_fx: 1.08 + dx0,
                lower_barrier_l: 1.02,
                upper_barrier_u: 1.15,
                domestic_rate: r,
                foreign_rate: 0.025 + dx2,
                atm_vol: 0.085 + dx3,
                risk_reversal_25d: -0.012,
                butterfly_25d: 0.0035,
                overhedge_survival_weight: 0.65,
                maturity_years: 0.5,
            },
        );
        let tmr = ql_exotics::evaluate_bouzoubaa_timer_option_variance_clock(
            &ql_exotics::BouzoubaaTimerOptionVarianceClockSpec {
                spot_price: spot,
                strike_price: 100.0,
                risk_free_rate: r,
                target_budget_vol: 0.16 + dx3,
                nominal_horizon_years: 1.0,
                market_implied_vol: vol,
                expected_realized_vol: 0.20,
                vol_of_vol_xi: 0.30,
            },
        );
        let bl = ql_portfolio::evaluate_meucci_black_litterman_woodbury_equivalence(
            &ql_portfolio::MeucciBlackLittermanWoodburySpec {
                covariance_3x3: [
                    [0.040 + dx0, 0.010, 0.004],
                    [0.010, 0.025 + dx1, 0.005],
                    [0.004, 0.005, 0.010 + dx2],
                ],
                market_cap_weights: [0.50, 0.30, 0.20],
                risk_aversion_delta: 2.5,
                tau_scalar: 0.05,
                view_pick_vector: [1.0, -1.0, 0.0],
                view_target_q: 0.045 + dx3,
                view_uncertainty_omega: 0.0015,
            },
        );
        let ac = ql_portfolio::evaluate_gueant_almgren_chriss_hyperbolic_liquidation(
            &ql_portfolio::GueantAlmgrenChrissHyperbolicSpec {
                initial_inventory_x0: 100.0 + dx0 * 10.0,
                horizon_years_t: 0.25,
                temporary_impact_eta: 0.002 + dx1 * 0.1,
                permanent_impact_gamma: 0.0005,
                asset_volatility_sigma: vol,
                risk_aversion_lambda: 0.15 + dx3,
            },
        );
        let enet = ql_dml::evaluate_hastie_zou_elastic_net_and_grouping(
            &ql_dml::HastieZouElasticNetGroupingSpec {
                collinear_rho01: 0.92,
                feature_target_corrs: [0.68 + dx0, 0.65 + dx1, 0.42 + dx2],
                lasso_penalty_lambda1: 0.08,
                ridge_penalty_lambda2: 0.25 + dx3,
            },
        );
        let nk = ql_models::evaluate_ljungqvist_new_keynesian_three_equation_ree(
            &ql_models::LjungqvistNewKeynesianReeSpec {
                risk_aversion_sigma: 1.5 + dx0,
                discount_factor_beta: 0.99,
                calvo_phillips_kappa: 0.15 + dx1,
                taylor_inflation_phi_pi: 1.50 + dx2,
                taylor_output_phi_y: 0.125 + dx3,
                natural_rate_persistence_rho: 0.80,
                natural_rate_shock_rn: 0.015 + credit_bump,
            },
        );
        nss.total_fisher_weil_duration
            + tips.floored_tips_fair_value * 0.01
            + rec.rfv_coupon_bond_fair_value * 0.01
            + slv.dupire_local_vol_delta
            + mva.central_clearing_xva_savings
            + ss.posterior_fitted_log_futures
            + dnt.smile_adjusted_dnt_pv
            + tmr.timer_call_fair_value * 0.1
            + bl.posterior_view_expected_return
            + ac.optimal_mean_variance_objective
            + enet.elastic_net_objective_value
            + nk.equilibrium_policy_rate
    }
    pub fn brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let hwt = ql_curves::evaluate_hagan_west_tension_convex_forward_curve(
            &ql_curves::HaganWestTensionCurveSpec {
                pillar_intervals: [1.0, 1.0, 2.0, 3.0],
                discrete_forwards: [r, r + 0.006, r + 0.012, r + 0.015],
                tension_lambda: 1.5,
                eval_coord_x: 0.35,
            },
        );
        let dl = ql_credit::evaluate_duffie_lando_noisy_accounting_credit_spread(
            &ql_credit::DuffieLandoNoisyCreditSpec {
                reported_log_distance_x_hat: 0.45 - credit_bump * 5.0,
                initial_log_distance_x0: 0.55,
                asset_volatility_sigma: vol,
                risk_free_rate: r,
                payout_rate_delta: 0.015,
                report_age_t0: 0.5,
                accounting_noise_std_a: 0.18 + credit_bump * 2.0,
                recovery_rate: 0.40,
                horizon_tau: 2.0,
            },
        );
        let pac = ql_fixed_income::evaluate_fabozzi_pac_companion_cmo_waterfall(
            &ql_fixed_income::FabozziPacCompanionCmoSpec {
                initial_pool_balance: 100.0,
                wac_coupon: 0.06,
                tranche_coupon: 0.055,
                psa_lower_band: 100.0,
                psa_upper_band: 250.0,
                realized_psa_speed: 165.0 + credit_bump * 500.0,
                discount_yield: r + 0.01,
            },
        );
        let fm = ql_portfolio::evaluate_campbell_fama_macbeth_and_shanken_eiv(
            &ql_portfolio::CampbellFamaMacbethShankenSpec {
                asset_mean_excess_returns: [0.045 + r * 0.1, 0.062, 0.078, 0.095 + inflation * 0.2],
                estimated_factor_betas: [0.70, 0.95, 1.20, 1.45],
                idiosyncratic_var_sigma_eps2: 0.015 + credit_bump * 0.1,
                factor_var_sigma_f2: 0.022,
                sample_periods_t: 120.0,
            },
        );
        let svi = ql_vol::evaluate_gatheral_jacquier_arbitrage_free_svi_density(
            &ql_vol::GatheralJacquierSviDensitySpec {
                forward_price: spot,
                strike_price: spot * 1.02,
                expiry_years: 1.0,
                svi_a: 0.015 + credit_bump * 0.1,
                svi_b: 0.12,
                svi_rho: -0.35,
                svi_m: 0.02,
                svi_sigma: vol,
            },
        );
        let bgk = ql_models::evaluate_broadie_glasserman_kou_discrete_barrier_lookback(
            &ql_models::BroadieGlassermanKouDiscreteBarrierSpec {
                spot_price: spot,
                strike_price: 100.0,
                down_barrier_h: 90.0,
                risk_free_rate: r,
                dividend_yield: 0.01,
                volatility: vol,
                maturity_years: 1.0,
                monitoring_steps_m: 50.0,
            },
        );
        let ci = ql_portfolio::evaluate_bouchaud_cross_impact_eigenmode_execution(
            &ql_portfolio::BouchaudCrossImpactEigenmodeSpec {
                parent_order_q1: 50.0,
                parent_order_q2: 30.0,
                self_impact_lambda1: 0.04 + credit_bump,
                self_impact_lambda2: 0.03,
                cross_impact_lambda12: 0.015,
                resilience_decay_rho: 0.60,
                risk_aversion_gamma: 0.10,
                asset_volatility: vol,
            },
        );
        let ks = ql_portfolio::evaluate_karnosky_singer_multicurrency_attribution(
            &ql_portfolio::KarnoskySingerMultiCurrencySpec {
                port_asset_weights: [0.45, 0.35, 0.20],
                bench_asset_weights: [0.40, 0.35, 0.25],
                port_fx_weights: [0.50, 0.30, 0.20],
                bench_fx_weights: [0.40, 0.35, 0.25],
                local_port_returns: [0.065 + r * 0.1, 0.048, 0.032],
                local_bench_returns: [0.058, 0.045, 0.035],
                local_eurodeposit_rates: [r, 0.028, 0.012],
                base_cash_rate: r,
                spot_fx_returns: [0.0, 0.015 - credit_bump, -0.012],
            },
        );
        let mjd = ql_levy::evaluate_cont_tankov_merton_jump_diffusion(
            &ql_levy::ContTankovMertonJumpDiffusionSpec {
                spot_price: spot,
                strike_price: 105.0,
                risk_free_rate: r,
                dividend_yield: 0.01,
                diffusion_vol_sigma: vol,
                jump_intensity_lambda: 0.75 + credit_bump * 5.0,
                log_jump_mean_mu: -0.06,
                log_jump_std_delta: 0.14,
                maturity_years: 1.0,
            },
        );
        let cv = ql_portfolio::evaluate_campbell_viceira_intertemporal_hedging_saa(
            &ql_portfolio::CampbellViceiraIntertemporalSaaSpec {
                current_expected_excess_return_xt: 0.055 + inflation * 0.2,
                long_run_excess_return_xbar: 0.050,
                predictor_persistence_phi: 0.88,
                return_volatility_sigma_r: vol,
                predictor_volatility_sigma_x: 0.035,
                return_predictor_corr_rho: -0.68,
                relative_risk_aversion_gamma: 4.0,
                wealth_consumption_rho_w: 0.96,
            },
        );
        let ohl = ql_models::evaluate_ohlson_penman_linear_information_dynamics(
            &ql_models::OhlsonPenmanLinearInfoDynamicsSpec {
                current_book_value_bt: 100.0,
                current_earnings_xt: 14.5 + inflation * 20.0,
                current_dividends_dt: 4.5,
                other_information_shock_nu: 1.8 - credit_bump * 10.0,
                cost_of_equity_re: 0.09 + r * 0.2,
                abnormal_earnings_persistence_omega: 0.62,
                other_info_persistence_gamma: 0.45,
            },
        );
        let rev = ql_ad::evaluate_griewank_revolve_checkpointed_bsb_pde_adjoint(
            &ql_ad::GriewankRevolveBsbAdjointSpec {
                spot_price: spot,
                strike_price: 100.0,
                risk_free_rate: r,
                sigma_min: 0.15,
                sigma_max: vol + 0.05,
                bsb_switching_steepness: 2.0,
                maturity_years: 0.5,
            },
        );

        hwt.zero_coupon_bond_price
            + dl.noisy_credit_spread_bps * 0.001
            + pac.pac_tranche_pv * 0.01
            + fm.factor_risk_premium_gamma1
            + svi.breeden_litzenberger_pdf
            + bgk.discrete_bgk_down_and_out_call * 0.01
            + ci.total_cross_impact_execution_cost * 0.01
            + ks.total_active_base_return
            + mjd.merton_jump_diffusion_call * 0.01
            + cv.total_strategic_equity_allocation
            + ohl.lid_intrinsic_equity_value * 0.01
            + rev.primal_bsb_option_value * 0.01
    }
    pub fn hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let xccy = ql_curves::evaluate_brigo_xccy_mtm_collateralized_basis_curve(
            &ql_curves::BrigoXccyMtmBasisCurveSpec {
                spot_fx_s0: 1.08,
                domestic_csa_rate_rd: r + 0.01,
                foreign_local_ois_rate_rf: r,
                fx_forward_quote_f0t: 1.11,
                mtm_xccy_par_basis_spread_bm: -0.0025 - credit_bump * 0.2,
                fx_vol_sigma_s: 0.09,
                foreign_rate_vol_sigma_rf: 0.008,
                fx_rate_corr_rho: -0.35,
                maturity_years_t: 2.0,
            },
        );
        let iopo = ql_fixed_income::evaluate_fabozzi_specified_pool_io_po_payup(
            &ql_fixed_income::FabozziSpecifiedPoolIoPoSpec {
                par_notional: 100.0,
                wac_coupon: 0.065,
                pass_through_coupon: 0.060,
                current_mortgage_rate: r + 0.015,
                discount_yield: r + 0.01,
                base_turnover_cpr: 0.06,
                max_refi_cpr: 0.42,
                s_curve_steepness_kappa: 120.0,
                spec_pool_burnout_dampener: 0.45,
            },
        );
        let mo = ql_credit::evaluate_lando_marshall_olkin_fatal_shock_basket(
            &ql_credit::LandoMarshallOlkinContagionSpec {
                idiosyncratic_intensity_lambda1: 0.018 + credit_bump,
                idiosyncratic_intensity_lambda2: 0.024 + credit_bump,
                fatal_common_shock_lambda12: 0.012 + credit_bump * 0.5,
                recovery_rate: 0.40,
                risk_free_rate: r,
                maturity_years_t: 5.0,
            },
        );
        let vix = ql_vol::evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(
            &ql_vol::GatheralBergomiVixVvixSpec {
                forward_variance_swap_var: vol * vol,
                vvix_vol_of_vol_omega: 0.85,
                vix_option_strike_k: vol * 0.95,
                risk_free_rate: r,
                expiry_years_t: 0.25,
            },
        );
        let wwr = ql_xva::evaluate_gregory_hull_white_wrong_way_risk_cva(
            &ql_xva::GregoryHullWhiteWwrCvaSpec {
                baseline_exposure_mu: 12.0,
                exposure_volatility_sigma_x: 18.0,
                baseline_hazard_rate_lambda0: 0.025 + credit_bump,
                log_intensity_vol_sigma_lambda: 0.65,
                exposure_credit_corr_rho: 0.55,
                recovery_rate: 0.40,
                risk_free_rate: r,
                maturity_years_t: 3.0,
            },
        );
        let gs = ql_commodities::evaluate_geman_gibson_schwartz_commodity_futures(
            &ql_commodities::GemanGibsonSchwartzSpec {
                spot_price_s0: spot * 0.8,
                convenience_yield_delta0: 0.14 + inflation,
                long_run_convenience_alpha: 0.08,
                mean_reversion_kappa: 1.4,
                spot_volatility_sigma1: vol + 0.1,
                convenience_volatility_sigma2: 0.26,
                spot_convenience_corr_rho: 0.65,
                risk_free_rate: r,
                maturity_years_tau: 1.5,
            },
        );
        let fva = ql_fx::evaluate_clark_fx_forward_volatility_agreement(
            &ql_fx::ClarkFxForwardVolAgreementSpec {
                spot_fx: 1.08,
                domestic_rate_rd: r + 0.01,
                foreign_rate_rf: r,
                front_expiry_t1: 0.5,
                back_expiry_t2: 1.5,
                atm_vol_t1: 0.085,
                atm_vol_t2: 0.098 + credit_bump,
                smile_convexity_adj_bps: -18.0,
            },
        );
        let ch = ql_exotics::evaluate_rubinstein_simple_chooser_option(
            &ql_exotics::RubinsteinChooserOptionSpec {
                spot_price: spot,
                strike_price: 100.0,
                risk_free_rate: r,
                dividend_yield: 0.02,
                volatility: vol,
                choice_date_tc: 0.4,
                final_maturity_t: 1.0,
            },
        );
        let flam = ql_portfolio::evaluate_grinold_kahn_flam_with_transfer_coefficient(
            &ql_portfolio::GrinoldKahnFlamTransferSpec {
                information_coefficient_ic: 0.06,
                nominal_breadth_n: 200.0,
                cross_sectional_signal_corr_rho: 0.12 + credit_bump,
                active_weights_std_adj: [0.45, 0.30, -0.20],
                unconstrained_optimal_weights_std: [0.55, 0.25, -0.35],
                active_risk_aversion_lambda: 2.5,
            },
        );
        let es = ql_portfolio::evaluate_meucci_expected_shortfall_euler_decomposition(
            &ql_portfolio::MeucciExpectedShortfallEulerSpec {
                weights: [0.45, 0.35, 0.20],
                expected_returns: [0.06, 0.05, 0.035],
                asset_volatilities: [vol, vol * 0.85, 0.10],
                pairwise_corr_rho: 0.35 + credit_bump,
                confidence_alpha: 0.975,
                normal_quantile_z_alpha: 1.959963984540054,
                portfolio_skewness_gamma1: -0.45,
                portfolio_excess_kurtosis_gamma2: 1.5,
            },
        );
        let js = ql_dml::evaluate_hastie_james_stein_empirical_bayes_sure(
            &ql_dml::HastieJamesSteinSureSpec {
                sample_means_y: [0.08 + inflation, 0.03, -0.02 - credit_bump, 0.11],
                grand_prior_target_mu0: 0.04,
                noise_variance_sigma2: 0.0025,
                singular_values_squared: [4.2, 2.5, 1.1, 0.4],
                ridge_penalty_lambda: 0.8,
            },
        );
        let cc = ql_models::evaluate_campbell_cochrane_external_habit_sdf(
            &ql_models::CampbellCochraneExternalHabitSpec {
                utility_curvature_gamma: 2.0,
                consumption_growth_mean_gc: 0.0189,
                consumption_growth_vol_sigmac: 0.015,
                log_surplus_persistence_phis: 0.87,
                subjective_discount_rate_rho_time: 0.0116,
                current_log_surplus_deviation_ds: -0.10 - credit_bump * 2.0,
            },
        );

        xccy.fx_implied_foreign_csa_df_pf
            + iopo.specified_pool_payup_32nds * 0.01
            + mo.second_to_default_par_spread_bps * 0.001
            + vix.black76_vix_call_price
            + wwr.wwr_adjusted_cva * 0.1
            + gs.gibson_schwartz_futures_price * 0.01
            + fva.forward_start_atm_straddle_pct
            + ch.rubinstein_chooser_price * 0.01
            + flam.expected_value_added_va
            + es.cornish_fisher_tail_es
            + js.stein_unbiased_risk_estimate_sure * 10.0
            + cc.conditional_equity_risk_premium
    }
    pub fn brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let ois = ql_curves::evaluate_brigo_ois_compounded_in_arrears_and_fallback(
            &ql_curves::BrigoOisInArrearsFallbackSpec {
                ois_df_start_p0t1: (-r).exp(),
                ois_df_end_p0t2: (-r * 1.5).exp(),
                accrual_year_fraction_delta: 0.5,
                short_rate_vol_sigma_r: 0.012 + credit_bump * 0.2,
                mean_reversion_a: 0.15,
                isda_historical_fallback_spread: 0.0026161,
                realized_accrued_rfr_rate: r + 0.005,
                elapsed_accrual_fraction_w: 0.4,
            },
        );
        let frn = ql_fixed_income::evaluate_fabozzi_collared_frn_discount_margin(
            &ql_fixed_income::FabozziCollaredFrnDiscountMarginSpec {
                par_face: 100.0,
                index_rate_r: r,
                quoted_margin_qm: 0.012,
                discount_margin_dm: 0.015 + credit_bump,
                cap_strike: 0.065,
                floor_strike: 0.035,
                rate_normal_vol_sigma: 0.011,
                time_to_next_reset_tau1: 0.25,
            },
        );
        let cdsw = ql_credit::evaluate_schonbucher_cds_swaption_and_fep(
            &ql_credit::SchonbucherCdsSwaptionFepSpec {
                hazard_rate_lambda: 0.022 + credit_bump,
                recovery_rate: 0.40,
                risk_free_rate: r,
                option_expiry_te: 1.0,
                swap_tenor_tm: 5.0,
                strike_spread_k: 0.013,
                spread_lognormal_vol_sigma: 0.45,
            },
        );
        let berg = ql_vol::evaluate_bergomi_two_factor_skew_term_structure(
            &ql_vol::BergomiTwoFactorSkewTermStructureSpec {
                flat_forward_variance_xi0: vol * vol,
                vol_of_vol_omega: 1.8,
                mixing_weight_theta: 0.35,
                fast_mean_reversion_k1: 6.5,
                slow_mean_reversion_k2: 0.35,
                factor_corr_rho12: 0.25,
                spot_fast_corr_rho_s1: -0.75,
                spot_slow_corr_rho_s2: -0.45,
                maturity_years_t: 3.0,
            },
        );
        let saccr = ql_xva::evaluate_gregory_sa_ccr_regulatory_ead(
            &ql_xva::GregorySaCcrRegulatoryEadSpec {
                trade_notional_d: 1000.0,
                current_mtm_value_v: -8.0 + credit_bump * 100.0,
                posted_collateral_c: 4.0,
                swap_start_years_s: 1.0,
                swap_end_years_e: 6.0,
                forward_rate_p: r + 0.005,
                option_strike_k: 0.035,
                supervisory_option_vol: 0.50,
                supervisory_factor_sf: 0.005,
            },
        );
        let spk = ql_commodities::evaluate_eydeland_clean_spark_spread_option(
            &ql_commodities::EydelandCleanSparkSpreadOptionSpec {
                forward_power_price_fp: spot * 0.68,
                forward_gas_price_fg: 7.2 + inflation * 10.0,
                plant_heat_rate_h: 7.5,
                fixed_vom_and_carbon_k: 8.5,
                power_volatility_sigma_p: vol + 0.2,
                gas_volatility_sigma_g: 0.36,
                power_gas_corr_rho: 0.72,
                risk_free_rate: r,
                expiry_years_t: 0.5,
            },
        );
        let prdc = ql_fx::evaluate_wystup_prdc_three_factor_coupon(
            &ql_fx::WystupPrdcThreeFactorCouponSpec {
                initial_spot_fx_s0: 105.0,
                current_spot_fx_st: 108.0 - credit_bump * 50.0,
                foreign_coupon_rate_cf: 0.055,
                domestic_strike_coupon_cd: 0.035,
                domestic_rate_rd: 0.010,
                foreign_rate_rf: r + 0.01,
                spot_fx_vol_sigma_s: 0.11,
                dom_rate_vol_sigma_rd: 0.007,
                for_rate_vol_sigma_rf: 0.012,
                corr_fx_dom_rho_sd: 0.25,
                corr_fx_for_rho_sf: -0.40,
                maturity_years_t: 5.0,
            },
        );
        let comp = ql_exotics::evaluate_geske_compound_call_and_put_on_call(
            &ql_exotics::GeskeCompoundOptionSpec {
                spot_price_s0: spot,
                compound_strike_k1: 6.5,
                underlying_strike_k2: 102.0,
                risk_free_rate: r,
                dividend_yield: 0.015,
                volatility: vol,
                compound_expiry_t1: 0.5,
                underlying_expiry_t2: 1.5,
            },
        );
        let mpir = ql_portfolio::evaluate_qian_multi_period_alpha_decay_ir(
            &ql_portfolio::QianMultiPeriodAlphaDecayIrSpec {
                initial_ic0: 0.075,
                alpha_decay_rate_kappa: 2.0 + credit_bump * 10.0,
                breadth_br: 180.0,
                proportional_tcost_rate_ctc: 0.08,
                active_risk_budget_sigma_a: 0.04,
            },
        );
        let enb = ql_portfolio::evaluate_roncalli_meucci_diversification_and_enb(
            &ql_portfolio::RoncalliMeucciDiversificationEnbSpec {
                weights: [0.30, 0.35, 0.35],
                asset_volatilities: [vol, vol * 0.8, 0.12],
                corr_01: 0.25 + credit_bump,
                corr_02: 0.10,
                corr_12: 0.30,
            },
        );
        let vecm = ql_dml::evaluate_tsay_johansen_vecm_cointegration(
            &ql_dml::TsayJohansenVecmCointegrationSpec {
                cointegrating_hedge_ratio_beta: 1.15,
                error_correction_speed_alpha1: -0.18 - credit_bump,
                error_correction_speed_alpha2: 0.08,
                current_spread_deviation_zt: 1.4,
                spread_innovation_vol_sigma_eps: 0.25,
                canonical_eigenvalues_lambda: [0.16, 0.03],
                sample_size_t: 250.0,
            },
        );
        let lrr = ql_models::evaluate_bansal_yaron_long_run_risks_sdf(
            &ql_models::BansalYaronLongRunRisksSpec {
                relative_risk_aversion_gamma: 7.5 + credit_bump * 20.0,
                elasticity_intertemporal_sub_psi: 1.5,
                campbell_shiller_kappa1: 0.997,
                growth_persistence_rho_x: 0.979,
                volatility_persistence_nu_sigma: 0.987,
                growth_shock_scale_phi_e: 0.044,
                baseline_consumption_vol_sigma: 0.027 + inflation * 0.1,
                vol_of_vol_sigma_w: 0.00023,
            },
        );

        ois.seasoned_period_blended_rfr_rate
            + frn.collared_frn_price * 0.01
            + cdsw.non_knockout_payer_swaption_bps * 0.001
            + berg.atmf_implied_vol_skew_slope.abs()
            + saccr.sa_ccr_exposure_at_default_ead * 0.01
            + spk.bjerksund_stensland_spark_call * 0.01
            + prdc.prdc_coupon_present_value
            + comp.call_on_call_price * 0.01
            + mpir.net_active_value_added
            + enb.shannon_entropy_effective_bets_enb * 0.1
            + vecm.mean_reversion_half_life_periods * 0.1
            + lrr.total_consumption_claim_risk_premium
    }
    pub fn brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let cms = ql_curves::evaluate_andersen_piterbarg_cms_replication_convexity(
            &ql_curves::AndersenPiterbargCmsReplicationSpec {
                forward_swap_rate_s0: r + 0.005,
                swap_tenor_years_n: 10.0,
                payment_frequency_m: 2.0,
                fixing_expiry_te: 5.0,
                normal_swaption_vol_sigma_n: 0.0085 + credit_bump * 0.1,
                ois_spread_delta_r: 0.0015,
                mean_reversion_speed_a: 0.04,
            },
        );
        let repo = ql_fixed_income::evaluate_tuckman_special_repo_fail_penalty_ctd(
            &ql_fixed_income::TuckmanSpecialRepoFailPenaltyCtdSpec {
                general_collateral_rate_gc: r + 0.01,
                tmpg_fail_penalty_rate_p: 0.030,
                issue_short_demand_ratio_d: 1.4 + credit_bump * 20.0,
                lendable_supply_elasticity_k: 0.55,
                bond_clean_price: 98.5,
                accrued_interest_at_spot: 0.85,
                coupon_rate_annual: 0.040,
                futures_invoice_price_at_delivery: 99.65,
                holding_period_year_frac_tau: 0.25,
            },
        );
        let cds_uf = ql_credit::evaluate_okane_isda_standard_cds_upfront_triangle(
            &ql_credit::OKaneIsdaStandardCdsUpfrontSpec {
                notional_amount: 1000.0,
                market_par_spread_s: 0.0185 + credit_bump,
                standard_running_coupon_c: 0.0100,
                recovery_rate_r: 0.40,
                ois_zero_rate_r: r,
                maturity_years_t: 5.0,
                accrued_fraction_alpha: 0.35,
            },
        );
        let ssr = ql_vol::evaluate_gatheral_ssr_tanaka_stop_loss_hedging(
            &ql_vol::GatheralSsrTanakaStopLossSpec {
                spot_price_s0: spot,
                strike_price_k: 100.0,
                risk_free_rate_r: r,
                atm_implied_vol_sigma0: vol,
                atm_skew_slope_dsigma_dk: -0.18,
                spot_vol_regression_beta: -0.24 - credit_bump * 2.0,
                maturity_years_t: 0.5,
                rebalancing_steps_n: 252.0,
            },
        );
        let kva = ql_xva::evaluate_gregory_saccr_kva_regulatory_capital(
            &ql_xva::GregorySaccrKvaRegulatoryCapitalSpec {
                initial_sa_ccr_ead0: 150.0 + credit_bump * 500.0,
                counterparty_risk_weight_rw: 0.75,
                pillar_1_capital_ratio_k: 0.105,
                cva_capital_add_on_multiplier: 0.35,
                shareholder_hurdle_rate_hk: 0.11,
                risk_free_discount_rate_r: r,
                exposure_amortization_decay_mu: 0.18,
                portfolio_maturity_years_t: 5.0,
            },
        );
        let crk = ql_commodities::evaluate_eydeland_crack_spread_three_asset_kirk(
            &ql_commodities::EydelandThreeAssetCrackSpreadSpec {
                crude_forward_barrel_fc: spot * 0.76,
                gasoline_forward_barrel_fg: spot * 0.94 + inflation * 20.0,
                distillate_forward_barrel_fh: spot * 0.98,
                fixed_refining_strike_k: 6.5,
                crude_vol_sigma_c: vol + 0.1,
                gasoline_vol_sigma_g: 0.38,
                distillate_vol_sigma_h: 0.36,
                corr_gas_dist_rho_gh: 0.82,
                corr_product_crude_rho_pc: 0.86,
                risk_free_rate_r: r,
                expiry_years_t: 0.5,
            },
        );
        let prdc_cap =
            ql_fx::evaluate_wystup_fx_prdc_coupon_convexity(&ql_fx::WystupFxPrdcCallableCapSpec {
                spot_fx_s0: 108.0 - credit_bump * 40.0,
                initial_fx_s_init: 105.0,
                foreign_coupon_cf: 0.055,
                domestic_deduction_cd: 0.035,
                cap_rate_k_cap: 0.030,
                domestic_rate_rd: 0.012,
                foreign_rate_rf: r + 0.01,
                fx_vol_sigma_s: 0.115,
                dom_rate_vol_sigma_d: 0.008,
                fx_dom_corr_rho_sd: 0.35,
                coupon_horizon_t: 7.0,
            });
        let eso = ql_exotics::evaluate_haug_eso_vesting_barrier_dilution(
            &ql_exotics::HaugEsoVestingBarrierDilutionSpec {
                spot_price_s0: spot,
                strike_price_k: 100.0,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                contractual_maturity_t: 10.0,
                vesting_period_tv: 3.0,
                forfeiture_hazard_lambda_f: 0.06 + credit_bump,
                suboptimal_exercise_multiple_m: 1.8,
                grant_to_shares_ratio_gamma: 0.05,
            },
        );
        let erc = ql_portfolio::evaluate_qian_risk_parity_spinu_log_barrier(
            &ql_portfolio::QianSpinuRiskParityLogBarrierSpec {
                asset_volatilities: [vol, vol * 0.7, 0.08 + credit_bump],
                risk_budgets_b: [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
                corr_01: 0.25,
                corr_02: -0.10,
                corr_12: 0.15,
            },
        );
        let cvar = ql_portfolio::evaluate_rockafellar_uryasev_cvar_minimization(
            &ql_portfolio::RockafellarUryasevCvarAuxiliarySpec {
                portfolio_weights: [0.6, 0.4],
                scenario_losses_asset0: [-0.02, 0.01, 0.04 + credit_bump, 0.09 + credit_bump],
                scenario_losses_asset1: [-0.01, 0.005, 0.02, 0.05],
                confidence_alpha: 0.75,
                smoothing_epsilon: 0.002,
            },
        );
        let pls =
            ql_dml::evaluate_htf_pls_vs_pcr_krylov_shrinkage(&ql_dml::HtfPlsKrylovShrinkageSpec {
                pc_eigenvalues_d2: [4.5 + credit_bump * 10.0, 1.8, 0.4],
                ols_pc_projections_u: [0.8, 1.4, 1.1],
                ridge_penalty_lambda: 0.6,
                noise_variance_sigma2: 0.05,
            });
        let fomc = ql_models::evaluate_piazzesi_swanson_fomc_futures_hike_tree(
            &ql_models::PiazzesiSwansonFomcFuturesHikeSpec {
                current_target_fed_funds_r0: r,
                futures_implied_month_avg_rate_f0: r + 0.0012 + inflation * 0.02,
                meeting_day_of_month_d: 15.0,
                total_days_in_month_n: 30.0,
                standard_hike_step_delta: 0.0025,
                piazzesi_swanson_term_premium_bps: 4.5 + credit_bump * 100.0,
                next_meeting_futures_rate_f1: r + 0.0038,
            },
        );

        cms.convexity_adjusted_cms_rate
            + repo.net_basis_after_specialness * 0.1
            + cds_uf.cash_settlement_upfront_amount * 0.01
            + ssr.minimum_variance_smile_delta
            + kva.lifetime_capital_valuation_adjustment_kva * 0.1
            + crk.crack_spread_call_price_per_bbl * 0.01
            + prdc_cap.capped_prdc_coupon_pv
            + eso.diluted_per_option_grant_cost * 0.01
            + erc.portfolio_volatility
            + cvar.minimized_auxiliary_cvar
            + pls.pls1_fitted_norm_sq * 0.05
            + fomc.physical_hike_probability_p
    }
    pub fn andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let ois_turn = ql_curves::evaluate_darbyshire_central_bank_step_ois_and_turn(
            &ql_curves::DarbyshireCentralBankStepOisTurnSpec {
                current_overnight_ois_r0: r,
                meeting_1_time_years_t1: 0.15,
                meeting_1_expected_step_bps: 22.0 + credit_bump * 200.0,
                meeting_2_time_years_t2: 0.35,
                meeting_2_expected_step_bps: 18.0,
                turn_of_year_spike_bps: 35.0,
                turn_effective_day_weight_tau: 3.0 / 360.0,
                horizon_years_t: 0.50,
            },
        );
        let cb_mw = ql_fixed_income::evaluate_despiegeleer_cb_soft_call_make_whole(
            &ql_fixed_income::DeSpiegeleerCbSoftCallMakeWholeSpec {
                spot_price_s0: spot,
                par_value: 100.0,
                conversion_ratio_cr: 1.0,
                soft_call_trigger_multiple: 1.30,
                make_whole_max_bonus_shares: 0.18,
                coupon_rate_c: 0.025,
                risky_bond_yield_y: r + 0.015 + credit_bump,
                risk_free_rate_r: r,
                dividend_and_borrow_cost_q: 0.018,
                equity_volatility_sigma: vol,
                maturity_years_t: 4.0,
            },
        );
        let cdx = ql_credit::evaluate_okane_cds_index_intrinsic_skew_basis(
            &ql_credit::OKaneCdsIndexIntrinsicSkewSpec {
                constituent_spreads_s: [
                    0.0060 + credit_bump * 0.5,
                    0.0120 + credit_bump,
                    0.0300 + credit_bump * 1.5,
                ],
                market_quoted_index_spread: 0.0155 + credit_bump,
                standard_index_coupon_c: 0.0100,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let vswp = ql_vol::evaluate_gatheral_varswap_log_strip_tail_truncation(
            &ql_vol::GatheralVarSwapTruncationJumpSpec {
                forward_price_f0: spot,
                diffusive_vol_sigma: vol,
                lower_cutoff_moneyness_kl_over_f: 0.65,
                upper_cutoff_moneyness_ku_over_f: 1.40,
                jump_intensity_lambda_j: 0.35 + credit_bump * 5.0,
                mean_jump_log_size_mu_j: -0.12,
                jump_log_vol_sigma_j: 0.08,
                maturity_years_t: 0.5,
            },
        );
        let fva_ov = ql_xva::evaluate_gregory_albanese_fva_fca_fba_overlap(
            &ql_xva::GregoryAlbaneseFvaOverlapSpec {
                expected_positive_exposure_epe: 42.0,
                expected_negative_exposure_ene: 28.0,
                bank_cds_hazard_rate_lambda_b: 0.015 + credit_bump,
                bank_recovery_rate_rb: 0.40,
                bond_cds_liquidity_basis_s_liq: 0.0045,
                treasury_deposit_spread_s_dep: 0.0035,
                counterparty_survival_weight: 0.94,
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let wkg = ql_commodities::evaluate_geman_working_curve_inventory_convenience_yield(
            &ql_commodities::GemanWorkingInventoryConvenienceSpec {
                spot_commodity_price_s0: spot * 0.82 + inflation * 20.0,
                normalized_inventory_i_over_icrit: 0.55,
                stockout_convenience_yield_y0: 0.32,
                abundance_convenience_floor_yinf: 0.01,
                working_curve_steepness_kappa: 1.8,
                physical_storage_cost_rate_c: 0.04,
                risk_free_rate_r: r,
                baseline_spot_vol_sigma0: vol,
                scarcity_vol_amplification_eta: 0.65,
                futures_maturity_t: 1.0,
            },
        );
        let fxfly = ql_fx::evaluate_wystup_malz_fx_broker_fly_smile_inversion(
            &ql_fx::WystupMalzFxSmileInversionSpec {
                spot_fx_s0: 1.085 - credit_bump * 0.5,
                domestic_rate_rd: r + 0.005,
                foreign_rate_rf: 0.030,
                atm_dns_vol_sigma_atm: 0.085 + credit_bump,
                risk_reversal_25d_rr: -0.012,
                butterfly_25d_bf: 0.0035,
                expiry_years_t: 0.5,
            },
        );
        let lkb = ql_exotics::evaluate_goldman_sosin_gatto_lookback_reflection(
            &ql_exotics::GoldmanSosinGattoLookbackSpec {
                spot_price_s0: spot,
                realized_minimum_smin: spot * 0.92,
                realized_maximum_smax: spot * 1.06,
                fixed_strike_k: 100.0,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let lw = ql_portfolio::evaluate_ledoit_wolf_constant_correlation_shrinkage(
            &ql_portfolio::LedoitWolfConstantCorrelationSpec {
                sample_volatilities: [vol, vol * 0.8, 0.14 + credit_bump],
                sample_correlations: [0.65, 0.15, 0.40],
                asymptotic_pi_hat: 0.045,
                asymptotic_rho_hat: 0.012,
                sample_observations_t: 60.0,
            },
        );
        let ac = ql_portfolio::evaluate_almgren_chriss_cara_optimal_execution(
            &ql_portfolio::AlmgrenChrissCaraExecutionSpec {
                initial_shares_x0: 100.0,
                liquidation_horizon_t: 1.0,
                asset_price_volatility_sigma: vol,
                temporary_impact_eta: 0.02 + credit_bump,
                permanent_impact_gamma: 0.005,
                cara_risk_aversion_lambda: 0.5,
            },
        );
        let estar = ql_dml::evaluate_tsay_star_threshold_autoregression_band(
            &ql_dml::TsayEstarThresholdAutoregressionSpec {
                current_mispricing_z: 1.85 + credit_bump * 10.0,
                transaction_cost_threshold_c: 0.90,
                transition_smoothness_gamma: 0.85,
                inner_band_reversion_speed_kappa_in: 0.04,
                outer_arbitrage_reversion_speed_kappa_out: 0.38,
                innovation_volatility_sigma_eps: 0.22,
            },
        );
        let acm = ql_models::evaluate_acm_cochrane_piazzesi_affine_term_premium(
            &ql_models::AcmCochranePiazzesiTermPremiumSpec {
                forward_rates_1y_to_5y: [
                    r,
                    r + 0.003,
                    r + 0.007 + inflation * 0.05,
                    r + 0.005,
                    r + 0.004,
                ],
                tent_weights_gamma: [-1.2, 0.4, 2.1, 0.5, -1.5],
                long_run_neutral_short_rate_r_star: 0.035,
                physical_var_mean_reversion_phi: 0.22,
                market_price_of_level_risk_lambda0: 0.35,
                market_price_of_slope_risk_lambda1: 0.80,
                horizon_years_n: 10.0,
            },
        );

        ois_turn.all_in_discount_factor_p0t
            + cb_mw.convertible_bond_fair_value * 0.01
            + cdx.rpv01_weighted_intrinsic_spread * 10.0
            + vswp.total_fair_variance_swap_strike
            + fva_ov.non_overlapping_net_fva_plus_dva * 0.1
            + wkg.storage_equilibrium_futures_price * 0.01
            + fxfly.inverted_25d_call_strike_k25c
            + lkb.floating_strike_lookback_call * 0.01
            + lw.equal_weight_shrunk_variance * 10.0
            + ac.cara_certainty_equivalent_cost * 0.01
            + estar.local_regime_half_life_periods * 0.1
            + acm.fitted_10y_nominal_yield * 10.0
    }
    pub fn darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let tb = ql_curves::evaluate_henrard_tenor_basis_multiplicative_spread(
            &ql_curves::HenrardTenorBasisMultiplicativeSpreadSpec {
                ois_discount_rate_r: r,
                forward_3m_period1_f1: r + 0.003,
                forward_3m_period2_f2: r + 0.005,
                quoted_forward_6m_f6: r + 0.0065 + credit_bump * 0.4,
                vol_3m_period1_sigma1: vol * 0.85,
                vol_3m_period2_sigma2: vol * 0.90,
                subperiod_correlation_rho12: 0.85,
                start_time_years_t0: 1.0,
            },
        );
        let msr = ql_fixed_income::evaluate_fabozzi_agency_msr_valuation(
            &ql_fixed_income::FabozziAgencyMsrValuationSpec {
                gross_servicing_fee_rate: 0.0025,
                ancillary_float_income_rate: 0.0004,
                annual_servicing_cost_rate: 0.0006,
                wac_mortgage_coupon: 0.060,
                current_primary_mortgage_rate: r + 0.015,
                base_turnover_cpr: 0.06,
                max_refi_cpr: 0.42,
                refi_scurve_steepness_kappa: 120.0,
                lender_recapture_fraction_alpha: 0.28,
                ois_discount_plus_oas_rate: r + 0.035 + credit_bump,
            },
        );
        let cdo = ql_credit::evaluate_okane_cdo_base_correlation_tranche(
            &ql_credit::OKaneCdoBaseCorrelationTrancheSpec {
                portfolio_hazard_rate_lambda: 0.018 + credit_bump * 1.2,
                recovery_rate_r: 0.40,
                attachment_point_k1: 0.03,
                detachment_point_k2: 0.07,
                base_correlation_k1_rho1: 0.24,
                base_correlation_k2_rho2: 0.36,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let fwd_sm = ql_vol::evaluate_gatheral_forward_start_cliquet_smile(
            &ql_vol::GatheralForwardStartCliquetSmileSpec {
                spot_price_s0: spot,
                reset_strike_moneyness_alpha: 1.02,
                reset_time_years_t1: 0.50,
                maturity_time_years_t2: 1.50,
                forward_initial_variance_v0: vol * vol,
                mean_reversion_speed_kappa: 1.8,
                vol_of_vol_eta: 0.55 + credit_bump * 5.0,
                spot_vol_correlation_rho: -0.65,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
            },
        );
        let csa_mta = ql_xva::evaluate_gregory_pykhtin_csa_threshold_mta_cva(
            &ql_xva::GregoryPykhtinCsaThresholdMtaCvaSpec {
                uncollateralized_exposure_vol_sigma_v: 35.0 + credit_bump * 100.0,
                current_mtm_value_v0: 12.0,
                csa_threshold_h: 15.0,
                minimum_transfer_amount_mta: 2.5,
                independent_amount_ia: 4.0,
                margin_period_of_risk_days: 10.0,
                counterparty_hazard_rate_lambda: 0.022 + credit_bump,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: r,
                horizon_years_t: 3.0,
            },
        );
        let swing = ql_commodities::evaluate_eydeland_gas_swing_take_or_pay(
            &ql_commodities::EydelandGasSwingTakeOrPaySpec {
                forward_gas_prices_f: [
                    3.20 + inflation * 5.0,
                    3.85 + inflation * 6.0,
                    3.45 + credit_bump * 8.0,
                ],
                contract_strike_price_k: 3.30,
                min_period_volume_qmin: 0.50,
                max_period_volume_qmax: 1.50,
                min_cumulative_acq_qmin_tot: 2.20,
                max_cumulative_acq_qmax_tot: 3.80,
                shortfall_makeup_penalty_per_unit: 0.85,
                spot_gas_volatility_sigma: 0.45,
                risk_free_rate_r: r,
            },
        );
        let kiko =
            ql_fx::evaluate_wystup_fx_window_kiko_barrier(&ql_fx::WystupFxWindowKikoBarrierSpec {
                spot_fx_s0: 1.085 - credit_bump * 0.5,
                strike_fx_k: 1.080,
                lower_knock_in_barrier_l: 1.020,
                upper_knock_out_barrier_u: 1.155,
                window_start_years_t1: 0.15,
                window_end_years_t2: 0.75,
                maturity_years_t: 1.00,
                domestic_rate_rd: r + 0.005,
                foreign_rate_rf: 0.030,
                volatility_sigma: 0.095 + credit_bump,
                upper_ko_rebate_r: 0.0080,
            });
        let stulz = ql_exotics::evaluate_stulz_two_asset_rainbow_option(
            &ql_exotics::StulzTwoAssetRainbowOptionSpec {
                spot_asset_1_s1: spot,
                spot_asset_2_s2: spot * 1.05,
                strike_price_k: 98.0,
                vol_asset_1_sigma1: vol,
                vol_asset_2_sigma2: vol * 1.15,
                correlation_rho: 0.45,
                dividend_yield_1_q1: 0.015,
                dividend_yield_2_q2: 0.020,
                risk_free_rate_r: r,
                maturity_years_t: 1.0,
            },
        );
        let jm = ql_portfolio::evaluate_jagannathan_ma_no_short_sale_shrinkage(
            &ql_portfolio::JagannathanMaCovarianceShrinkageSpec {
                asset_volatilities: [vol * 1.2, vol * 0.9, 0.14 + credit_bump],
                correlation_12: 0.85,
                correlation_13: 0.25,
                correlation_23: 0.30,
                weight_upper_cap_u: 0.65,
            },
        );
        let ow = ql_portfolio::evaluate_obizhaeva_wang_transient_lob_execution(
            &ql_portfolio::ObizhaevaWangTransientLobExecutionSpec {
                total_order_size_x0: 100.0,
                execution_horizon_t: 1.0,
                lob_depth_density_q0: 25.0,
                permanent_impact_gamma: 0.008,
                lob_resilience_speed_rho: 3.2 - credit_bump * 10.0,
                asset_volatility_sigma: vol,
                risk_aversion_lambda: 0.15,
            },
        );
        let hac = ql_dml::evaluate_newey_west_andrews_hac_variance(
            &ql_dml::NeweyWestAndrewsHacVarianceSpec {
                sample_variance_gamma0: 0.040 + credit_bump,
                autocovariance_lag1_gamma1: 0.022,
                autocovariance_lag2_gamma2: 0.012,
                autocovariance_lag3_gamma3: 0.006,
                sample_size_n: 250.0,
                ols_point_estimate_beta: 0.035,
            },
        );
        let gab = ql_models::evaluate_gabaix_rare_disaster_resilience(
            &ql_models::GabaixRareDisasterResilienceSpec {
                disaster_probability_p: 0.018 + credit_bump * 0.5,
                consumption_recovery_in_disaster_b: 0.72,
                dividend_recovery_in_disaster_f: 0.55,
                crra_risk_aversion_gamma: 4.0,
                time_preference_rate_rho: 0.025,
                normal_growth_rate_g: 0.020,
                resilience_deviation_h_hat: -0.008 - credit_bump,
                resilience_mean_reversion_phi_h: 0.20,
            },
        );

        tb.ois_discounted_basis_swaplet_pv * 0.1
            + msr.recapture_adjusted_msr_fair_value_pct
            + cdo.mezzanine_fair_running_spread_bps * 0.01
            + fwd_sm.stochastic_vol_forward_start_call * 0.1
            + csa_mta.collateralized_cva
            + swing.total_swing_contract_value
            + kiko.window_kiko_call_value * 10.0
            + stulz.call_on_max_value * 0.1
            + jm.constrained_gmv_portfolio_volatility * 10.0
            + ow.risk_adjusted_total_cost * 0.01
            + hac.hac_robust_t_statistic
            + gab.linearity_generating_price_dividend_ratio * 0.1
    }
    pub fn henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let dswap = ql_curves::evaluate_brigo_differential_swap_quanto_fra(
            &ql_curves::BrigoDifferentialSwapQuantoFraSpec {
                domestic_ois_rate_rd: r,
                domestic_forward_rate_ld: r + 0.003,
                foreign_forward_rate_lf: 0.032 + credit_bump * 0.5,
                foreign_rate_vol_sigma_lf: vol,
                fx_spot_vol_sigma_fx: 0.12,
                rate_fx_correlation_rho: -0.45,
                reset_time_years_t1: 2.0,
                accrual_fraction_tau: 0.50,
            },
        );
        let muni = ql_fixed_income::evaluate_fabozzi_callable_muni_ytw_refunding(
            &ql_fixed_income::FabozziCallableMuniYtwRefundingSpec {
                tax_exempt_coupon_rate_c: 0.050,
                muni_mmd_benchmark_yield_ym: r * 0.78 + credit_bump * 0.5,
                treasury_yield_y_ust: r + 0.008,
                statutory_marginal_tax_rate_tau: 0.37,
                call_protection_years_tc: 5.0,
                final_maturity_years_t: 10.0,
                call_price_par: 100.0,
                muni_yield_volatility_sigma: 0.18 + credit_bump,
            },
        );
        let step_cb = ql_credit::evaluate_bielecki_rating_step_up_coupon_bond(
            &ql_credit::BieleckiRatingStepUpCouponBondSpec {
                base_coupon_rate_c0: 0.045,
                downgrade_step_up_spread_dc: 0.0125,
                ig_to_hy_downgrade_intensity_lambda_down: 0.08 + credit_bump * 2.0,
                hy_to_ig_upgrade_intensity_lambda_up: 0.03,
                ig_default_hazard_rate_lambda_d_ig: 0.005 + credit_bump * 0.5,
                hy_default_hazard_rate_lambda_d_hy: 0.035 + credit_bump,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let bkm = ql_vol::evaluate_bakshi_kapadia_madan_moment_strip(
            &ql_vol::BakshiKapadiaMadanMomentStripSpec {
                spot_price_s0: spot,
                atm_implied_vol_sigma0: vol,
                smile_skew_slope_beta: -0.25 - credit_bump * 2.0,
                smile_curvature_gamma: 0.40,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                maturity_years_t: 0.25,
            },
        );
        let brk_cva = ql_xva::evaluate_gregory_break_clause_extinguishable_cva(
            &ql_xva::GregoryBreakClauseExtinguishableCvaSpec {
                unmitigated_swap_epe_profile: [18.0, 32.0, 45.0],
                break_clause_date_years_tb: 2.0,
                final_maturity_years_t: 5.0,
                counterparty_hazard_rate_lambda: 0.025 + credit_bump,
                break_execution_probability_p_exec: 0.85,
                break_settlement_slippage_fraction: 0.08,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: r,
            },
        );
        let pwr = ql_commodities::evaluate_bessembinder_lemmon_power_forward_premium(
            &ql_commodities::BessembinderLemmonPowerForwardPremiumSpec {
                baseline_marginal_fuel_cost_c0: 28.0 + inflation * 50.0,
                expected_load_to_capacity_ratio_u: 0.84,
                load_shock_volatility_sigma_u: 0.12 + credit_bump,
                supply_stack_convexity_exponent_gamma: 4.5,
                peaking_merit_slope_alpha: 65.0,
                producer_variance_aversion_a_var: 0.015,
                retailer_spike_skew_preference_b_skew: 0.0045,
            },
        );
        let fader =
            ql_fx::evaluate_wystup_fx_fader_range_accrual(&ql_fx::WystupFxFaderRangeAccrualSpec {
                spot_fx_s0: 1.085 - credit_bump * 0.5,
                strike_fx_k: 1.080,
                accrual_corridor_lower_l: 1.040,
                accrual_corridor_upper_u: 1.130,
                domestic_rate_rd: r + 0.005,
                foreign_rate_rf: 0.030,
                volatility_sigma: 0.095 + credit_bump,
                maturity_years_t: 1.0,
            });
        let ext_opt = ql_exotics::evaluate_longstaff_extendible_option(
            &ql_exotics::LongstaffExtendibleOptionSpec {
                spot_price_s0: spot,
                initial_strike_k1: 102.0,
                extended_strike_k2: 105.0,
                initial_maturity_t1: 0.50,
                extended_maturity_t2: 1.25,
                holder_extension_fee_a: 1.50,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
            },
        );
        let kz = ql_portfolio::evaluate_kan_zhou_three_fund_combine(
            &ql_portfolio::KanZhouThreeFundCombineSpec {
                sample_sharpe_squared_theta2: 0.28 + credit_bump,
                gmv_sharpe_squared_theta_g2: 0.12,
                equal_weight_sharpe_squared_theta_ew2: 0.16,
                number_of_assets_n: 25.0,
                sample_observations_t: 120.0,
                risk_aversion_gamma: 3.0,
            },
        );
        let rbsa = ql_portfolio::evaluate_sharpe_return_based_style_analysis(
            &ql_portfolio::SharpeReturnBasedStyleAnalysisSpec {
                fund_mean_return_mu_p: 0.105 - credit_bump,
                fund_return_variance_var_p: 0.032,
                style_benchmark_means: [0.095, 0.075, 0.045],
                style_benchmark_variances: [0.040, 0.028, 0.010],
                fund_style_covariances: [0.024, 0.011, -0.002],
            },
        );
        let dm_hln = ql_dml::evaluate_diebold_mariano_hln_forecast_test(
            &ql_dml::DieboldMarianoHlnForecastTestSpec {
                model1_mse_loss_l1: 0.025 + credit_bump * 0.2,
                model2_mse_loss_l2: 0.036,
                forecast_error_correlation_rho12: 0.45,
                loss_differential_autocov_lag1: 0.0004,
                forecast_horizon_steps_h: 3.0,
                evaluation_sample_size_t: 80.0,
            },
        );
        let lw_rstar = ql_models::evaluate_laubach_williams_natural_rate_r_star(
            &ql_models::LaubachWilliamsNaturalRateRStarSpec {
                prior_trend_growth_g0: 0.018,
                prior_other_factor_z0: -0.004 - credit_bump * 0.2,
                nominal_policy_rate_i: r + 0.015,
                expected_inflation_pi_e: inflation,
                realized_core_inflation_pi: inflation + 0.004,
                is_curve_surprise_dy: 0.006,
                is_rate_sensitivity_a_r: 0.09,
                phillips_slope_b_y: 0.15,
                signal_to_noise_ratio_lambda_z: 0.042,
            },
        );

        dswap.domestic_discounted_diff_swaplet_pv * 0.1
            + muni.callable_muni_market_price * 0.01
            + step_cb.step_up_covenanted_bond_value * 0.01
            + bkm.cboe_skew_index_level * 0.01
            + brk_cva.break_clause_mitigated_cva
            + pwr.equilibrium_power_forward_price * 0.01
            + fader.fade_in_call_value * 10.0
            + ext_opt.writer_extendible_call_value * 0.1
            + kz.three_fund_expected_out_of_sample_utility * 10.0
            + rbsa.active_information_ratio
            + dm_hln.hln_corrected_dm_statistic * 0.1
            + lw_rstar.taylor_rule_implied_nominal_rate * 10.0
    }
    pub fn brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let sofr_fut = ql_curves::evaluate_andersen_sofr_futures_fra_convexity(
            &ql_curves::AndersenSofrFuturesFraConvexitySpec {
                unadjusted_fra_forward_rate_f0: r + 0.005,
                short_rate_volatility_sigma_r: 0.0115 + credit_bump * 0.2,
                hull_white_mean_reversion_kappa: 0.06,
                accrual_start_years_t1: 3.0,
                accrual_fraction_tau: 0.25,
                contract_notional_mm: 10.0,
            },
        );
        let mbs_prep = ql_fixed_income::evaluate_hayre_richard_roll_mbs_prepayment(
            &ql_fixed_income::HayreRichardRollMbsPrepaymentSpec {
                pool_wac_coupon_rate: 0.062,
                prevailing_mortgage_rate_r_mtg: r + 0.015,
                loan_age_months_wala: 24.0,
                cumulative_burnout_exposure: 0.35,
                burnout_sensitivity_beta: 0.60,
                s_curve_steepness_gamma: 65.0,
                seasonal_turnover_multiplier: 1.08,
            },
        );
        let cdx_opt = ql_credit::evaluate_pedersen_okane_index_option_fep_smile(
            &ql_credit::PedersenOkaneIndexOptionFepSmileSpec {
                spot_index_spread_s0: 0.0085 + credit_bump,
                strike_spread_k: 0.0090,
                standard_running_coupon_c: 0.0100,
                spread_volatility_sigma: 0.45 + credit_bump * 2.0,
                option_expiry_years_te: 0.50,
                index_maturity_years_tm: 5.50,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: r,
            },
        );
        let dup_lv = ql_vol::evaluate_gatheral_total_variance_dupire_local_vol(
            &ql_vol::GatheralTotalVarianceDupireLocalVolSpec {
                log_forward_moneyness_y: -0.08,
                atm_implied_vol_sigma0: vol,
                log_moneyness_skew_s: -0.24 - credit_bump,
                log_moneyness_convexity_c: 0.35,
                term_structure_slope_theta_t: 0.015,
                maturity_years_t: 0.50,
            },
        );
        let simm_mva = ql_xva::evaluate_gregory_kenyon_simm_mva_funding(
            &ql_xva::GregoryKenyonSimmMvaFundingSpec {
                initial_simm_margin_im0: 42.0,
                ccp_cleared_margin_ratio: 0.62,
                bank_funding_spread_s_f: 0.0095 + credit_bump * 0.5,
                custodian_remuneration_spread_s_c: 0.0015,
                bank_hazard_rate_lambda_b: 0.012,
                counterparty_hazard_rate_lambda_c: 0.018 + credit_bump,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let rss_stor = ql_commodities::evaluate_routledge_seppi_spatt_stockout_backwardation(
            &ql_commodities::RoutledgeSeppiSpattStockoutBackwardationSpec {
                spot_commodity_price_s0: 78.0 + inflation * 100.0,
                discretionary_inventory_i: 0.24,
                critical_stockout_threshold_i_star: 0.20,
                base_abundant_volatility_sigma0: 0.22 + credit_bump,
                stockout_volatility_spike_eta: 0.85,
                interest_plus_storage_cost_rate_c: r + 0.025,
                convenience_scarcity_slope_alpha: 0.16,
                maturity_years_t: 0.50,
            },
        );
        let cip = ql_fx::evaluate_baba_packer_cip_cross_currency_basis_wedge(
            &ql_fx::BabaPackerCipCrossCurrencyBasisWedgeSpec {
                spot_fx_s0: 1.0800,
                market_outright_forward_fx_f: 1.1015 - credit_bump * 0.5,
                domestic_usd_ois_rate_rd: r + 0.01,
                foreign_ois_rate_rf: 0.0310,
                gsib_leverage_ratio_capital_charge_bps: 14.0,
                fx_swap_bid_ask_and_funding_haircut_bps: 4.5,
                maturity_years_t: 1.0,
            },
        );
        let pcs_bar = ql_exotics::evaluate_carr_chou_put_call_symmetry_barrier_replication(
            &ql_exotics::CarrChouPutCallSymmetryBarrierReplicationSpec {
                spot_price_s0: spot,
                strike_price_k: 102.0,
                down_and_out_barrier_h: 90.0,
                risk_free_rate_r: r,
                dividend_yield_q: 0.04,
                volatility_sigma: vol,
                maturity_years_t: 0.75,
            },
        );
        let ps_shrink = ql_portfolio::evaluate_pastor_stambaugh_factor_prior_shrinkage(
            &ql_portfolio::PastorStambaughFactorPriorShrinkageSpec {
                sample_ols_alpha_hat: 0.035 - credit_bump * 0.5,
                sample_factor_beta_hat: 1.12,
                residual_return_variance_sigma_eps2: 0.018 + credit_bump,
                prior_mispricing_std_sigma_alpha: 0.020,
                benchmark_excess_mean_mu_m: 0.065,
                benchmark_return_variance_sigma_m2: 0.025,
                sample_length_periods_t: 60.0,
                risk_aversion_gamma: 3.0,
            },
        );
        let tb_app = ql_portfolio::evaluate_treynor_black_appraisal_portfolio(
            &ql_portfolio::TreynorBlackAppraisalPortfolioSpec {
                security_alphas: [0.024, 0.018, 0.012],
                security_market_betas: [1.15, 0.92, 1.05],
                security_residual_variances: [0.035 + credit_bump, 0.025, 0.020],
                market_excess_return_mu_m: 0.065,
                market_return_variance_sigma_m2: 0.0256,
            },
        );
        let eg_ecm = ql_dml::evaluate_engle_granger_cointegration_ecm_half_life(
            &ql_dml::EngleGrangerCointegrationEcmHalfLifeSpec {
                variance_asset_x: 0.042,
                covariance_xy: 0.038,
                spread_lag1_autocorrelation_phi: 0.86 + credit_bump,
                spread_innovation_sigma_e: 0.015,
                current_disequilibrium_spread_z: 0.025,
                sample_observations_t: 250.0,
                time_step_years_dt: 1.0 / 252.0,
            },
        );
        let ez_sdf = ql_models::evaluate_epstein_zin_weil_recursive_utility_sdf(
            &ql_models::EpsteinZinWeilRecursiveUtilitySdfSpec {
                subjective_discount_beta: 0.988,
                relative_risk_aversion_gamma: 7.5 + credit_bump * 10.0,
                intertemporal_elasticity_psi: 1.5,
                expected_consumption_growth_gc: 0.018 - inflation * 0.1,
                consumption_growth_volatility_sigma_c: 0.022,
                wealth_portfolio_return_volatility_sigma_w: 0.145,
                consumption_wealth_correlation_rho_cw: 0.65,
            },
        );

        sofr_fut.hull_white_futures_fra_convexity_bps
            + mbs_prep.annualized_cpr * 10.0
            + cdx_opt.payer_index_swaption_pv_bps * 0.1
            + dup_lv.dupire_local_volatility * 10.0
            + simm_mva.bilateral_simm_mva_charge
            + rss_stor.equilibrium_forward_price * 0.1
            + cip.covered_interest_parity_basis_bps * 0.1
            + pcs_bar.static_replicated_down_and_out_call * 0.1
            + ps_shrink.optimal_active_satellite_weight * 10.0
            + tb_app.combined_tangency_sharpe_ratio
            + eg_ecm.engle_granger_adf_t_statistic.abs() * 0.1
            + ez_sdf.hansen_jagannathan_max_sharpe_bound * 10.0
    }
    pub fn andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric(&self) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let pca_fly = ql_curves::evaluate_litterman_scheinkman_butterfly_pca_carry_rolldown(
            &ql_curves::LittermanScheinkmanButterflyPcaSpec {
                short_wing_duration_d1: 1.92,
                belly_duration_d2: 4.55,
                long_wing_duration_d3: 8.40,
                short_wing_slope_loading_s1: -0.65,
                belly_slope_loading_s2: -0.10,
                long_wing_slope_loading_s3: 0.48,
                short_wing_convexity_c1: 4.2,
                belly_convexity_c2: 23.5,
                long_wing_convexity_c3: 82.0,
                parallel_yield_shift_bps: 50.0 + credit_bump * 500.0,
            },
        );
        let tf_cb = ql_fixed_income::evaluate_tsiveriotis_fernandes_split_discount_convertible(
            &ql_fixed_income::TsiveriotisFernandesSplitConvertibleSpec {
                stock_price_s0: spot,
                conversion_ratio_cr: 1.0,
                par_redemption_value: 100.0,
                annual_coupon_rate_c: 0.035,
                risk_free_rate_r: r,
                credit_spread_sc: 0.025 + credit_bump,
                dividend_yield_q: 0.012,
                equity_volatility_sigma: vol,
                maturity_years_t: 3.0,
            },
        );
        let dt_cdo = ql_credit::evaluate_hull_white_double_t_cdo_tail_copula(
            &ql_credit::HullWhiteDoubleTCdoTailCopulaSpec {
                portfolio_default_probability_q: 0.045 + credit_bump,
                recovery_rate_r: 0.40,
                factor_correlation_rho: 0.30,
                student_t_degrees_of_freedom_nu: 4.0,
                equity_detachment_k1: 0.03,
                senior_attachment_k2: 0.15,
                maturity_years_t: 5.0,
            },
        );
        let lee_wing = ql_vol::evaluate_roger_lee_wing_moment_slope_bound(
            &ql_vol::RogerLeeWingMomentSlopeBoundSpec {
                critical_upper_moment_p_tilde: 3.5,
                critical_lower_moment_q_tilde: 1.8,
                atm_total_variance_w0: vol * vol * 0.5,
                extreme_log_moneyness_k: -0.45 - credit_bump,
                maturity_years_t: 0.50,
            },
        );
        let cds_cont = ql_xva::evaluate_brigo_capponi_cds_contagion_gap_cva(
            &ql_xva::BrigoCapponiCdsContagionGapCvaSpec {
                reference_hazard_rate_lambda_ref: 0.025 + credit_bump,
                protection_seller_hazard_rate_lambda_c: 0.020 + credit_bump * 0.5,
                contractual_running_spread_s_cds: 0.010,
                contagion_spread_jump_gamma: 0.85,
                simultaneous_default_prob_p_sim: 0.06,
                reference_recovery_rate_r_ref: 0.40,
                counterparty_recovery_rate_r_c: 0.40,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let sch_fut = ql_commodities::evaluate_schwartz_one_factor_mean_reverting_log_futures(
            &ql_commodities::SchwartzOneFactorLogFuturesSpec {
                spot_commodity_price_s0: 76.0 + inflation * 100.0,
                risk_neutral_long_run_log_price_mu_star: 4.25,
                mean_reversion_speed_kappa: 0.85,
                spot_log_volatility_sigma: 0.34 + credit_bump,
                option_expiry_years_te: 0.50,
                futures_maturity_years_tf: 1.00,
                option_strike_price_k: 74.0,
                risk_free_rate_r: r,
            },
        );
        let dcd = ql_fx::evaluate_wystup_dual_currency_deposit_yield_enhancement(
            &ql_fx::WystupDualCurrencyDepositYieldSpec {
                spot_fx_s0: 1.0850 - credit_bump * 0.5,
                conversion_strike_fx_k: 1.0700,
                domestic_deposit_rate_rd: r + 0.005,
                foreign_deposit_rate_rf: 0.030,
                fx_volatility_sigma: 0.105 + credit_bump,
                deposit_tenor_years_t: 0.25,
                bank_structuring_margin_pct: 0.15,
            },
        );
        let pass_opt = ql_exotics::evaluate_hyer_lipton_pugachevsky_passport_option(
            &ql_exotics::HyerLiptonPugachevskyPassportOptionSpec {
                spot_price_s0: spot,
                initial_trading_account_w0: 2.0,
                max_position_limit_q_max: 1.0,
                risk_free_rate_r: r,
                cost_of_carry_b: 0.02,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let grs = ql_portfolio::evaluate_gibbons_ross_shanken_multivariate_factor_ftest(
            &ql_portfolio::GibbonsRossShankenFactorFTestSpec {
                alpha_quadratic_form_sr_alpha2: 0.145 + credit_bump,
                benchmark_factor_sharpe_ratio_sr_f: 0.42,
                sample_observations_t: 180.0,
                test_assets_count_n: 25.0,
                factor_count_k: 3.0,
            },
        );
        let sqrt_imp = ql_portfolio::evaluate_bouchaud_gatheral_square_root_metaorder_impact(
            &ql_portfolio::BouchaudGatheralSquareRootImpactSpec {
                metaorder_shares_q: 250_000.0,
                average_daily_volume_adv: 5_000_000.0,
                daily_return_volatility_sigma_d: 0.018 + credit_bump * 0.2,
                y_prefactor_constant: 0.65,
                impact_concavity_exponent_delta: 0.50,
                propagator_decay_exponent_beta: 0.55,
                post_trade_elapsed_ratio_tau_over_t: 1.0,
            },
        );
        let bns_har = ql_dml::evaluate_barndorff_nielsen_shephard_bipower_jump_har_cj(
            &ql_dml::BarndorffNielsenShephardHarCjSpec {
                daily_realized_variance_rv: 0.00042 + credit_bump * 1e-3,
                daily_bipower_variation_bv: 0.00028,
                daily_tripower_quarticity_tq: 1.15e-7,
                weekly_continuous_bv_w: 0.00024,
                monthly_continuous_bv_m: 0.00020,
                intraday_sampling_frequency_m: 78.0,
            },
        );
        let ll_cay = ql_models::evaluate_lettau_ludvigson_cay_consumption_wealth_forecast(
            &ql_models::LettauLudvigsonCayForecastSpec {
                log_consumption_c: 9.865,
                log_asset_wealth_a: 11.20 - credit_bump,
                log_labor_income_y: 9.25,
                asset_wealth_cointegrating_share_omega: 0.31,
                unconditional_equity_premium_gamma0: 0.062,
                cay_predictive_slope_beta_cay: 2.15,
                conditional_market_volatility_sigma_m: 0.16 + credit_bump,
            },
        );

        pca_fly.convexity_pnl_advantage_bps
            + tf_cb.tsiveriotis_fernandes_convertible_value * 0.1
            + dt_cdo.senior_tranche_tail_spread_bps * 0.1
            + lee_wing.asymptotic_wing_implied_volatility * 10.0
            + cds_cont.total_contagion_adjusted_cds_cva * 10.0
            + sch_fut.futures_call_option_price * 0.1
            + dcd.enhanced_dcd_coupon_rate * 10.0
            + pass_opt.carry_adjusted_passport_value * 0.1
            + grs.grs_f_statistic
            + sqrt_imp.post_trade_decayed_impact_bps * 0.1
            + bns_har.har_cj_annualized_volatility_forecast * 10.0
            + ll_cay.conditional_market_sharpe_ratio * 10.0
    }
    pub fn litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let pit_opt = ql_curves::evaluate_piterbarg_collateral_choice_option_funding_wedge(
            &ql_curves::PiterbargCollateralOptionFundingWedgeSpec {
                spot_price_s0: spot,
                strike_price_k: 100.0,
                csa_collateral_rate_r_c: r,
                equity_repo_financing_rate_r_r: r + 0.004,
                unsecured_funding_spread_s_f: 0.015 + credit_bump,
                dividend_yield_q: 0.015,
                equity_volatility_sigma_s: vol,
                funding_spread_volatility_sigma_f: 0.008,
                equity_funding_correlation_rho_sf: -0.35,
                maturity_years_t: 2.0,
            },
        );
        let tuck_roll = ql_fixed_income::evaluate_tuckman_par_spot_forward_rolldown_horizon_return(
            &ql_fixed_income::TuckmanRolldownHorizonReturnSpec {
                initial_maturity_years_t: 7.0,
                holding_horizon_years_dt: 1.0,
                current_bond_yield_y_t: r + 0.012,
                rolled_maturity_yield_y_t_minus_dt: r + 0.009,
                short_financing_rate_r_dt: r + 0.004,
                annual_coupon_rate_c: 0.045,
                modified_duration_at_horizon_d_mod: 5.20,
                convexity_at_horizon_conv: 32.5,
                expected_parallel_shift_bps: 15.0 + credit_bump * 200.0,
            },
        );
        let kmv_dd = ql_credit::evaluate_merton_kmv_distance_to_default_credit_equity_hedge(
            &ql_credit::MertonKmvDistanceToDefaultHedgeSpec {
                asset_value_v_a: 165.0 - credit_bump * 500.0,
                short_term_debt_std: 45.0,
                long_term_debt_ltd: 90.0,
                physical_asset_drift_mu_a: 0.075,
                risk_free_rate_r: r,
                asset_volatility_sigma_a: vol,
                horizon_years_t: 3.0,
            },
        );
        let bart_sabr = ql_vol::evaluate_hagan_bartlett_normal_sabr_minimum_variance_greeks(
            &ql_vol::HaganBartlettNormalSabrGreeksSpec {
                forward_rate_f: r + 0.008,
                strike_rate_k: 0.042,
                sabr_alpha: 0.045 + credit_bump,
                sabr_beta: 0.50,
                sabr_rho: -0.35,
                sabr_nu: 0.40,
                option_expiry_t: 2.0,
                ois_discount_factor_df: (-2.0 * r).exp(),
            },
        );
        let mpor_spk = ql_xva::evaluate_pykhtin_brownian_bridge_mpor_cashflow_spike_epe(
            &ql_xva::PykhtinMporCashflowSpikeEpeSpec {
                annualized_portfolio_mtm_volatility_sigma: 28.0 + credit_bump * 100.0,
                margin_period_of_risk_years_tau_mpor: 10.0 / 252.0,
                csa_threshold_plus_mta: 2.5,
                independent_amount_held_ia: 1.0,
                contractual_cashflow_paid_in_mpor: 6.0,
                post_payment_lookback_weight_theta: 0.60,
                counterparty_hazard_rate_lambda: 0.022 + credit_bump,
                counterparty_recovery_rate_r: 0.40,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let ftr_opt =
            ql_commodities::evaluate_deng_johnson_locational_transmission_congestion_option(
                &ql_commodities::DengJohnsonTransmissionCongestionOptionSpec {
                    sink_hub_forward_price_fb: 68.0 + inflation * 100.0,
                    source_hub_forward_price_fa: 54.0,
                    transmission_loss_factor_eta: 0.96,
                    wheeling_tariff_strike_k_trans: 4.5,
                    sink_hub_volatility_sigma_b: 0.48 + credit_bump,
                    source_hub_volatility_sigma_a: 0.42,
                    hub_price_correlation_rho_ab: 0.78,
                    risk_free_rate_r: r,
                    expiry_years_t: 0.50,
                },
            );
        let fx_dig = ql_fx::evaluate_wystup_fx_digital_call_spread_skew_overhedge(
            &ql_fx::WystupFxDigitalSkewOverhedgeSpec {
                spot_fx_s0: 1.0850 - credit_bump * 0.5,
                digital_strike_k: 1.1000,
                domestic_rate_rd: r + 0.005,
                foreign_rate_rf: 0.030,
                atm_or_strike_vol_sigma_k: 0.095 + credit_bump,
                smile_slope_dsigma_dk: -0.18,
                call_spread_half_width_dk: 0.0050,
                expiry_years_t: 0.50,
            },
        );
        let cv_lb = ql_exotics::evaluate_conze_viswanathan_fixed_strike_lookback_excursion(
            &ql_exotics::ConzeViswanathanFixedStrikeLookbackSpec {
                spot_price_s0: spot,
                fixed_strike_k: 105.0,
                risk_free_rate_r: r + 0.02,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let hot_es = ql_portfolio::evaluate_tasche_litterman_hot_spots_best_hedge_es_decomposition(
            &ql_portfolio::TascheLittermanHotSpotsBestHedgeSpec {
                focal_position_weight_w1: 0.45,
                core_portfolio_weight_w2: 0.55,
                focal_asset_volatility_sigma1: vol + credit_bump,
                core_portfolio_volatility_sigma2: 0.16,
                asset_core_correlation_rho12: 0.50,
                focal_expected_return_mu1: 0.08,
                core_expected_return_mu2: 0.06,
                gaussian_es_tail_multiplier_m_alpha: 2.338,
            },
        );
        let bab = ql_portfolio::evaluate_frazzini_pedersen_betting_against_beta_zero_beta_capm(
            &ql_portfolio::FrazziniPedersenBettingAgainstBetaSpec {
                low_beta_portfolio_beta_l: 0.65,
                high_beta_portfolio_beta_h: 1.45,
                low_beta_expected_return_r_l: 0.088,
                high_beta_expected_return_r_h: 0.118,
                risk_free_rate_r_f: r,
                shadow_funding_constraint_spread_psi: 0.022 + credit_bump,
                bab_factor_volatility_sigma_bab: 0.105,
            },
        );
        let pt_dir = ql_dml::evaluate_pesaran_timmermann_directional_timing_hit_ratio_test(
            &ql_dml::PesaranTimmermannDirectionalTimingSpec {
                true_positive_up_fraction_p11: 0.34 - credit_bump,
                true_negative_down_fraction_p00: 0.25,
                false_positive_up_fraction_p01: 0.20,
                sample_observations_t: 250.0,
                mean_absolute_return_move_bps: 85.0 + credit_bump * 500.0,
            },
        );
        let lmp = ql_models::evaluate_lucas_mehra_prescott_endowment_tree_equity_premium(
            &ql_models::LucasMehraPrescottEndowmentTreeSpec {
                mean_consumption_growth_mu_c: 0.018 - inflation * 0.1,
                consumption_growth_vol_sigma_c: 0.036 + credit_bump,
                dividend_consumption_leverage_lambda: 3.0,
                crra_risk_aversion_gamma: 2.5 + credit_bump * 10.0,
                subjective_discount_factor_beta: 0.99,
                historical_observed_equity_premium: 0.0618,
            },
        );

        pit_opt.total_piterbarg_funding_adjusted_call_price * 0.1
            + tuck_roll.horizon_expected_excess_return_pct
            + (kmv_dd.structural_credit_spread_bps * 0.01 + kmv_dd.physical_distance_to_default_dd)
            + bart_sabr.bartlett_minimum_variance_delta * 10.0
            + mpor_spk.total_mpor_collateralized_cva
            + ftr_opt.ftr_option_value_per_mwh * 0.1
            + fx_dig.skew_adjusted_digital_call_price * 10.0
            + cv_lb.fixed_strike_lookback_call_value * 0.1
            + hot_es.total_portfolio_expected_shortfall * 10.0
            + (bab.bab_factor_sharpe_ratio + bab.theoretical_funding_constrained_bab_alpha)
            + pt_dir.pesaran_timmermann_z_statistic
            + lmp.lucas_tree_price_dividend_ratio * 0.1
    }
    pub fn piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let gsw = ql_curves::evaluate_gurkaynak_sack_wright_svensson_treasury_forward_curve(
            &ql_curves::GswSvenssonTreasuryCurveSpec {
                svensson_beta0: r + 0.012,
                svensson_beta1: -0.012,
                svensson_beta2: 0.018 + credit_bump,
                svensson_beta3: -0.010,
                svensson_tau1: 1.80,
                svensson_tau2: 6.50,
                maturity_years_t: 10.0,
                on_the_run_par_yield_y_otr: r + 0.008,
            },
        );
        let mbs_esf = ql_fixed_income::evaluate_fabozzi_hayre_mbs_excess_servicing_gfee_strip(
            &ql_fixed_income::FabozziHayreMbsExcessServicingGfeeSpec {
                pool_balance_b0: 100.0,
                gross_wac_rate_w: 0.0640,
                net_passthrough_coupon_c: 0.0550,
                agency_guaranty_fee_g: 0.0045 + credit_bump * 0.1,
                minimum_servicing_fee_s_min: 0.0025,
                prepayment_speed_cpr: 0.12 - credit_bump * 2.0,
                scheduled_principal_amortization_a: 0.018,
                discount_rate_r: r + 0.01,
                pool_horizon_years_t: 10.0,
            },
        );
        let irb_ga = ql_credit::evaluate_vasicek_gordy_asrf_irb_capital_granularity_adjustment(
            &ql_credit::VasicekGordyAsrfIrbGranularitySpec {
                unconditional_pd: 0.015 + credit_bump,
                downturn_lgd: 0.45,
                asset_correlation_rho: 0.18,
                confidence_level_alpha: 0.999,
                effective_maturity_years_m: 3.0,
                herfindahl_concentration_index_h: 0.025,
            },
        );
        let cw_vol = ql_vol::evaluate_carr_wu_proportional_volatility_smile_drift(
            &ql_vol::CarrWuProportionalVolSmileDriftSpec {
                forward_price_f: spot,
                strike_price_k: 95.0,
                atm_volatility_sigma0: vol,
                implied_vol_of_vol_omega: 0.65 + credit_bump * 5.0,
                spot_vol_correlation_rho: -0.70,
                implied_vol_drift_mu: -0.04,
                time_to_expiry_tau: 0.50,
            },
        );
        let greg_cva = ql_xva::evaluate_gregory_incremental_marginal_netting_set_cva(
            &ql_xva::GregoryIncrementalMarginalNettingCvaSpec {
                existing_portfolio_mtm_m1: 12.0,
                new_trade_mtm_m2: -4.0,
                existing_portfolio_vol_sigma1: 18.0 + credit_bump * 50.0,
                new_trade_vol_sigma2: 10.0,
                exposure_correlation_rho12: -0.35,
                counterparty_hazard_rate_lambda: 0.024 + credit_bump,
                counterparty_recovery_rate_r: 0.40,
                ois_discount_rate_r: r,
                maturity_years_t: 4.0,
            },
        );
        let wb_stor = ql_commodities::evaluate_working_brennan_storage_convenience_yield_curve(
            &ql_commodities::WorkingBrennanStorageConvenienceYieldSpec {
                spot_price_s0: 75.0 + inflation * 100.0,
                inventory_to_consumption_ratio_i: 0.11,
                minimum_pipeline_inventory_i_min: 0.08,
                convenience_yield_convexity_alpha: 0.018,
                stockout_steepness_exponent_gamma: 0.85,
                physical_storage_cost_rate_u: 0.045,
                risk_free_rate_r: r,
                futures_maturity_years_t: 1.0,
            },
        );
        let fx_fe = ql_fx::evaluate_wystup_fx_forward_extra_participating_collar(
            &ql_fx::WystupFxForwardExtraCollarSpec {
                spot_fx_s0: 1.0850 - credit_bump * 0.5,
                worst_case_strike_k: 1.0950,
                downside_trigger_barrier_b: 1.0350,
                domestic_rate_rd: r + 0.005,
                foreign_rate_rf: 0.030,
                strike_vol_sigma_k: 0.095 + credit_bump,
                put_wing_skew_slope_xi: -0.14,
                expiry_years_t: 0.50,
            },
        );
        let ex_on_ex = ql_exotics::evaluate_rubinstein_margrabe_exchange_on_exchange_option(
            &ql_exotics::RubinsteinMargrabeExchangeOnExchangeSpec {
                spot_asset1_s1: spot * 1.05,
                spot_asset2_s2: 100.0,
                compound_strike_fraction_q: 0.04,
                dividend_yield_q1: 0.015,
                dividend_yield_q2: 0.025,
                volatility_sigma1: vol + 0.02,
                volatility_sigma2: vol,
                asset_correlation_rho: 0.55,
                underlying_expiry_years_t2: 1.50,
            },
        );
        let hlz_fdr = ql_portfolio::evaluate_harvey_liu_zhu_factor_tstat_fdr_haircut(
            &ql_portfolio::HarveyLiuZhuFactorFdrHaircutSpec {
                observed_factor_t_stat: 3.25 - credit_bump * 10.0,
                sample_months_t: 240.0,
                total_tested_factors_m: 316.0,
                target_fdr_alpha: 0.05,
                prior_true_factor_prob_p0: 0.15,
                prior_true_sharpe_std_sigma0: 0.35,
            },
        );
        let hs_kyle = ql_portfolio::evaluate_holden_subrahmanyam_multi_insider_kyle_equilibrium(
            &ql_portfolio::HoldenSubrahmanyamMultiInsiderKyleSpec {
                number_of_insiders_n: 4.0,
                fundamental_value_vol_sigma_v: 2.50 + credit_bump * 10.0,
                noise_trader_volume_vol_sigma_u: 120.0,
                private_valuation_surprise_dv: 1.80,
                prior_reference_price_p0: spot,
            },
        );
        let garch_ts = ql_dml::evaluate_bollerslev_tsay_garch_term_structure_kurtosis(
            &ql_dml::BollerslevTsayGarchTermStructureSpec {
                garch_intercept_omega: 0.000_004,
                arch_shock_alpha: 0.085,
                garch_persistence_beta: 0.895,
                current_daily_variance_h1: 0.000_32 + credit_bump * 0.005,
                innovation_kurtosis_kappa_z: 3.50,
                forecast_horizon_days_k: 21.0,
            },
        );
        let cs_cape = ql_models::evaluate_campbell_shiller_cape_excess_yield_forecast(
            &ql_models::CampbellShillerCapeExcessYieldSpec {
                real_equity_index_price_p: 4800.0 - credit_bump * 5000.0,
                ten_year_real_average_earnings_e10: 160.0,
                historical_mean_log_cape: 3.05,
                annual_valuation_reversion_speed_phi: 0.065,
                expected_real_earnings_growth_g: 0.018,
                dividend_payout_ratio_d: 0.45,
                ten_year_real_bond_yield_y10: r - inflation,
            },
        );

        gsw.off_the_run_liquidity_spread_bps * 0.1
            + mbs_esf.excess_servicing_strip_pv
            + irb_ga.total_granularity_adjusted_capital * 10.0
            + cw_vol.carr_wu_arbitrage_free_implied_vol * 10.0
            + greg_cva.incremental_trade_cva
            + wb_stor.equilibrium_futures_price * 0.1
            + fx_fe.forward_extra_client_npv * 10.0
            + ex_on_ex.exchange_on_exchange_call_value * 0.1
            + hlz_fdr.posterior_haircutted_sharpe_ratio
            + (hs_kyle.expected_clearing_price + hs_kyle.aggregate_insider_expected_profit * 0.01)
                * 0.1
            + garch_ts.annualized_term_structure_volatility * 10.0
            + cs_cape.implied_ten_year_equity_risk_premium * 100.0
    }
    pub fn gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let kw_shadow = ql_curves::evaluate_kim_wright_shadow_rate_survey_term_premium(
            &ql_curves::KimWrightShadowRateTermPremiumSpec {
                current_shadow_short_rate_s0: r - 0.015,
                effective_lower_bound_r_min: 0.00125,
                mean_reversion_speed_kappa: 0.22,
                risk_neutral_long_run_mean_theta_q: r + 0.014,
                physical_long_run_mean_theta_p: r + 0.002,
                shadow_rate_volatility_sigma_s: 0.014 + credit_bump * 0.2,
                survey_long_run_short_rate_r_surv: 0.028,
                survey_kalman_shrinkage_weight_w: 0.35,
                horizon_years_tau: 10.0,
            },
        );
        let bs_cb = ql_fixed_income::evaluate_brennan_schwartz_callable_puttable_convertible_game(
            &ql_fixed_income::BrennanSchwartzCallablePuttableCbSpec {
                stock_price_s0: spot,
                conversion_ratio_n_c: 1.0,
                par_notional_n: 100.0,
                coupon_rate_c: 0.025,
                risk_free_rate_r: r,
                credit_spread_s_c: 0.022 + credit_bump,
                dividend_yield_q: 0.012,
                equity_volatility_sigma_s: vol + 0.10,
                investor_put_price_k_put: 96.0,
                issuer_call_trigger_price_k_call: 130.0,
                maturity_years_t: 5.0,
            },
        );
        let ds_rmv = ql_credit::evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis(
            &ql_credit::DuffieSingletonRmvRfvBasisSpec {
                risk_free_rate_r: r,
                hazard_rate_lambda: 0.03 + credit_bump,
                recovery_rate_rec: 0.40,
                bond_coupon_rate_c: 0.055,
                deliverable_recovery_dispersion_sigma_r: 0.08 + credit_bump,
                maturity_years_t: 5.0,
            },
        );
        let fn_gamma = ql_vol::evaluate_fukasawa_neuberger_gamma_skewness_swap_strip(
            &ql_vol::FukasawaNeubergerGammaSkewnessSwapSpec {
                forward_price_f0: spot,
                atm_volatility_sigma0: vol,
                log_moneyness_skew_slope_alpha: -0.18 - credit_bump,
                log_moneyness_curvature_beta: 0.25,
                realized_spot_weighted_variance: 0.042,
                maturity_years_t: 0.50,
                variance_notional_n_var: 100.0,
            },
        );
        let bk_xva = ql_xva::evaluate_burgard_kjaer_semi_replication_fva_cva_dva(
            &ql_xva::BurgardKjaerSemiReplicationXvaSpec {
                expected_positive_exposure_epe: 25.0 + credit_bump * 50.0,
                expected_negative_exposure_ene: 14.0,
                counterparty_hazard_rate_lambda_c: 0.022 + credit_bump,
                counterparty_lgd_c: 0.60,
                bank_hazard_rate_lambda_b: 0.012,
                bank_lgd_b: 0.60,
                bank_funding_shortfall_fraction_gamma_f: 0.65,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let cd_spark = ql_commodities::evaluate_carmona_durrleman_spark_spread_lower_bound(
            &ql_commodities::CarmonaDurrlemanSparkSpreadSpec {
                electricity_forward_f_e: 62.0 + inflation * 100.0,
                natural_gas_forward_f_g: 7.20,
                heat_rate_mmbtu_per_mwh_h: 7.50,
                variable_om_strike_k_om: 6.50,
                electricity_volatility_sigma_e: 0.48 + credit_bump,
                gas_volatility_sigma_g: 0.42,
                power_gas_correlation_rho: 0.65,
                risk_free_rate_r: r,
                maturity_years_t: 0.75,
            },
        );
        let fx_tarf = ql_fx::evaluate_wystup_lipton_fx_tarf_target_redemption_cap(
            &ql_fx::WystupLiptonFxTarfSpec {
                spot_fx_s0: 1.09 - credit_bump * 0.5,
                tarf_strike_k: 1.06,
                target_redemption_cap_t_cap: 0.25,
                downside_gearing_leverage_l: 2.0,
                total_fixing_dates_m: 12.0,
                domestic_rate_rd: r,
                foreign_rate_rf: 0.030,
                fx_volatility_sigma: 0.10 + credit_bump,
                maturity_years_t: 1.0,
            },
        );
        let bgk_lb = ql_exotics::evaluate_broadie_glasserman_kou_discrete_lookback_option(
            &ql_exotics::BroadieGlassermanKouDiscreteLookbackSpec {
                spot_price_s0: spot,
                observed_running_minimum_m0: spot * 0.92,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                discrete_monitoring_dates_n: 52.0,
                maturity_years_t: 1.0,
            },
        );
        let nm_ff5 = ql_portfolio::evaluate_novy_marx_fama_french_five_factor_spanning_alpha(
            &ql_portfolio::NovyMarxFamaFrenchSpanningSpec {
                gross_profitability_premium_rmw: 0.042 + credit_bump * 0.2,
                value_premium_hml: 0.038,
                investment_premium_cma: 0.032,
                hml_on_rmw_loading_phi_rmw: 0.22,
                hml_on_cma_loading_phi_cma: 0.55,
                hml_spanning_residual_vol_sigma_e: 0.068,
                rmw_factor_vol_sigma_rmw: 0.085,
                hml_factor_vol_sigma_hml: 0.105,
                rmw_hml_correlation_rho: -0.42,
            },
        );
        let acerbi_srm = ql_portfolio::evaluate_acerbi_tasche_exponential_spectral_risk_measure(
            &ql_portfolio::AcerbiTascheExponentialSpectralRiskSpec {
                portfolio_expected_return_mu: 0.08 - credit_bump,
                portfolio_volatility_sigma: vol,
                return_skewness_gamma1: -0.55,
                return_excess_kurtosis_gamma2: 2.40,
                cara_spectral_risk_aversion_k: 25.0,
                reference_es_confidence_alpha: 0.975,
            },
        );
        let corsi_har = ql_dml::evaluate_corsi_har_realized_volatility_multi_horizon_forecast(
            &ql_dml::CorsiHarRealizedVolatilitySpec {
                har_intercept_c: 0.000_012,
                daily_rv_weight_beta_d: 0.36,
                weekly_rv_weight_beta_w: 0.34,
                monthly_rv_weight_beta_m: 0.25,
                current_daily_rv: 0.000_42 + credit_bump * 0.005,
                current_weekly_rv: 0.000_31,
                current_monthly_rv: 0.000_26,
                realized_quarticity_rq: 0.000_000_45,
                forecast_horizon_days_h: 22.0,
            },
        );
        let bll_nn = ql_models::evaluate_blanchard_lhuillier_lorenzoni_news_noise_kalman_filter(
            &ql_models::BlanchardLhuillierLorenzoniNewsNoiseSpec {
                fundamental_productivity_vol_sigma_a: 0.015,
                exogenous_news_noise_vol_sigma_e: 0.020 + credit_bump * 0.2,
                shannon_channel_capacity_bits_kappa: 1.50,
                observed_macro_news_signal_s: 0.025 + inflation * 0.2,
                consumption_demand_elasticity_gamma_c: 0.85,
                phillips_curve_slope_kappa_pc: 0.18,
            },
        );

        kw_shadow.survey_anchored_forward_term_premium * 100.0
            + kw_shadow.zlb_call_option_wedge_bps * 0.1
            + bs_cb.game_convertible_fair_value * 0.1
            + ds_rmv.rfv_risky_coupon_bond_price
            + ds_rmv.mod_r_restructuring_cds_basis_bps * 0.1
            + fn_gamma.gamma_swap_mark_to_market_npv
            + bk_xva.semi_replication_total_xva_charge
            + cd_spark.carmona_durrleman_spark_call_value
            + fx_tarf.client_tarf_net_present_value * 10.0
            + bgk_lb.bgk_discrete_lookback_call_value * 0.1
            + nm_ff5.combined_value_profitability_sharpe
            + acerbi_srm.exponential_spectral_risk_measure * 10.0
            + corsi_har.cumulative_average_h_day_rv * 1000.0
            + bll_nn.phillips_curve_inflation_impulse * 100.0
    }
    pub fn kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let jps_macro =
            ql_curves::evaluate_joslin_priebsch_singleton_unspanned_macro_term_structure(
                &ql_curves::JoslinPriebschSingletonUnspannedMacroSpec {
                    yield_curve_level_l: r + 0.012,
                    yield_curve_slope_s: -0.014,
                    yield_curve_curvature_c: 0.008,
                    real_output_growth_g: 0.024 - credit_bump,
                    core_inflation_gap_pi: inflation,
                    macro_on_yields_projection_beta: 0.25,
                    spanned_slope_risk_price_lambda_p: 0.18,
                    unspanned_macro_risk_price_lambda_m: -0.22,
                    q_measure_decay_kappa_q: 0.45,
                    bond_maturity_years_n: 10.0,
                },
            );
        let kf_ref = ql_fixed_income::evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency(
            &ql_fixed_income::KalotayFabozziRefundingEfficiencySpec {
                notional_par_n: 100.0,
                existing_coupon_rate_c_old: 0.055,
                refunding_market_yield_y_ref: r + 0.008,
                call_strike_price_k_call: 101.0,
                refunding_issuance_cost_c_iss: 0.65,
                yield_volatility_sigma_y: 0.22 + credit_bump,
                deferral_option_horizon_t_c: 1.50,
                remaining_maturity_years_t: 10.0,
            },
        );
        let hw_ntd = ql_credit::evaluate_hull_white_nth_to_default_basket_cds_spread(
            &ql_credit::HullWhiteNthToDefaultBasketSpec {
                basket_size_m: 5.0,
                individual_hazard_rate_lambda: 0.02 + credit_bump,
                recovery_rate_rec: 0.40,
                gaussian_copula_correlation_rho: 0.30 + credit_bump,
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let derman_reg = ql_vol::evaluate_derman_sticky_smile_regime_delta_adjustment(
            &ql_vol::DermanStickySmileRegimeSpec {
                forward_price_f: spot,
                strike_price_k: 95.0,
                atm_volatility_sigma0: vol,
                log_moneyness_skew_slope_b1: -0.15 - credit_bump,
                log_moneyness_curvature_b2: 0.20,
                risk_free_rate_r: r,
                maturity_years_t: 0.50,
            },
        );
        let hw_wwr = ql_xva::evaluate_hull_white_exponential_hazard_wrong_way_risk_cva(
            &ql_xva::HullWhiteExponentialWwrCvaSpec {
                portfolio_expected_mtm_mu_v: 8.0,
                portfolio_exposure_vol_sigma_v: 12.0 + credit_bump * 50.0,
                baseline_cds_hazard_rate_lambda0: 0.025 + credit_bump,
                wwr_exponential_sensitivity_b: 0.035,
                counterparty_lgd: 0.60,
                ois_discount_rate_r: r,
                maturity_years_t: 4.0,
            },
        );
        let gs_sam = ql_commodities::evaluate_gibson_schwartz_two_factor_futures_samuelson_vol(
            &ql_commodities::GibsonSchwartzSamuelsonFuturesSpec {
                spot_price_s0: 78.0 + inflation * 100.0,
                initial_convenience_yield_delta0: 0.14,
                convenience_yield_mean_reversion_kappa: 1.35,
                risk_neutral_long_run_convenience_alpha: 0.08,
                spot_volatility_sigma_s: 0.36 + credit_bump,
                convenience_yield_volatility_sigma_d: 0.28,
                spot_convenience_correlation_rho: 0.80,
                risk_free_rate_r: r,
                maturity_years_tau: 1.50,
            },
        );
        let quanto_fx = ql_fx::evaluate_reiner_wystup_quanto_fx_correlation_convexity(
            &ql_fx::ReinerWystupQuantoFxSpec {
                foreign_asset_spot_s0: spot,
                strike_price_k: 102.0,
                fixed_quanto_exchange_rate_q_fx: 1.0,
                domestic_rate_rd: r + 0.005,
                foreign_rate_rf: 0.025,
                foreign_dividend_yield_q: 0.015,
                asset_volatility_sigma_s: vol,
                fx_volatility_sigma_x: 0.11 + credit_bump,
                asset_fx_correlation_rho_sx: 0.45,
                maturity_years_t: 1.0,
            },
        );
        let tw_asian = ql_exotics::evaluate_turnbull_wakeman_edgeworth_arithmetic_asian_option(
            &ql_exotics::TurnbullWakemanArithmeticAsianSpec {
                spot_price_s0: spot,
                strike_price_k: 100.0,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let amp_vme = ql_portfolio::evaluate_asness_moskowitz_pedersen_value_momentum_everywhere(
            &ql_portfolio::AsnessMoskowitzPedersenVmeSpec {
                value_expected_excess_return_mu_v: 0.045 - credit_bump,
                momentum_expected_excess_return_mu_m: 0.062,
                value_volatility_sigma_v: 0.10 + credit_bump,
                momentum_volatility_sigma_m: 0.12,
                value_momentum_correlation_rho_vm: -0.52,
                tsmom_lookback_autocovariance_gamma: 0.18,
            },
        );
        let gp_aim = ql_portfolio::evaluate_garleanu_pedersen_dynamic_aim_portfolio_execution(
            &ql_portfolio::GarleanuPedersenAimPortfolioSpec {
                current_position_x_prev: 0.40,
                expected_alpha_return_f: 0.035 - credit_bump,
                alpha_mean_reversion_speed_phi: 0.25,
                return_variance_sigma2: 0.04,
                risk_aversion_gamma: 2.0,
                proportional_tcost_lambda: 1.50 + credit_bump * 10.0,
                discount_rate_rho_d: r,
            },
        );
        let bn_icp2 = ql_dml::evaluate_bai_ng_approximate_factor_icp2_diffusion_index(
            &ql_dml::BaiNgApproximateFactorIcp2Spec {
                cross_section_dimension_n: 120.0,
                time_series_dimension_t: 240.0,
                candidate_factor_count_k: 4.0,
                total_panel_variance_sigma_x2: 1.0 + credit_bump * 5.0,
                explained_variance_fraction_r2_k: 0.58,
                diffusion_forecast_beta_f: 0.35,
                factor_shock_signal_s_f: 1.40,
            },
        );
        let bgg_acc = ql_models::evaluate_bernanke_gertler_gilchrist_financial_accelerator(
            &ql_models::BernankeGertlerGilchristAcceleratorSpec {
                capital_asset_value_qk: 200.0,
                entrepreneur_net_worth_n: 80.0 - credit_bump * 500.0,
                csv_auditing_cost_fraction_mu: 0.12,
                efp_leverage_elasticity_nu: 0.08,
                risk_free_real_rate_r: r - inflation,
                unexpected_asset_price_shock_dq: 0.03,
            },
        );

        jps_macro.total_expected_bond_excess_return * 100.0
            + kf_ref.kalotay_refunding_efficiency_ratio * 10.0
            + hw_ntd.first_to_default_fair_spread_bps * 0.01
            + derman_reg.sticky_moneyness_delta * 10.0
            + hw_wwr.wrong_way_risk_adjusted_cva
            + gs_sam.gibson_schwartz_futures_price * 0.1
            + quanto_fx.quanto_call_domestic_value * 0.1
            + tw_asian.turnbull_wakeman_arithmetic_call_value * 0.1
            + amp_vme.combined_vme_sharpe_ratio
            + gp_aim.optimal_next_period_holding_x * 10.0
            + bn_icp2.bai_ng_pc_p2_criterion
            + bgg_acc.required_expected_return_on_capital * 100.0
    }
    pub fn jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let lm_cap = ql_curves::evaluate_lyashenko_mercurio_fmm_in_arrears_caplet_convexity(
            &ql_curves::LyashenkoMercurioFmmInArrearsCapletSpec {
                notional_n: 1_000.0,
                forward_rfr_rate_f0: r + 0.012,
                caplet_strike_k: 0.040,
                fixing_start_years_t0: 0.50,
                accrual_fraction_tau: 0.25,
                forward_vol_sigma_f: 0.24 + credit_bump,
                ois_discount_rate_r: r,
            },
        );
        let cpi_seas = ql_fixed_income::evaluate_tuckman_cpi_inflation_swap_seasonality_carry(
            &ql_fixed_income::TuckmanCpiInflationSwapSeasonalitySpec {
                notional_n: 100.0,
                nominal_treasury_yield_y_nom: r + 0.013,
                tips_real_yield_r_real: 0.019,
                inflation_risk_premium_irp: 0.0025 + credit_bump * 0.2,
                tips_liquidity_discount_liq: 0.0015,
                monthly_seasonal_cpi_factor_s_m: 0.0020 + inflation * 0.05,
                swap_maturity_years_t: 5.0,
            },
        );
        let ba_cir = ql_credit::evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds(
            &ql_credit::BrigoAlfonsiShiftedCirStochasticIntensitySpec {
                initial_intensity_y0: 0.022 + credit_bump,
                mean_reversion_speed_kappa: 0.45,
                long_run_intensity_theta: 0.028 + credit_bump,
                intensity_volatility_sigma: 0.14,
                deterministic_shift_psi: 0.004,
                rates_credit_correlation_rho: -0.35,
                short_rate_volatility_sigma_r: 0.012,
                risk_free_rate_r: r,
                recovery_rate_rec: 0.40,
                maturity_years_t: 5.0,
            },
        );
        let fpss_sv = ql_vol::evaluate_fouque_papanicolaou_sircar_multiscale_sv_smile(
            &ql_vol::FouquePapanicolaouSircarMultiscaleSvSpec {
                spot_price_s0: spot,
                strike_price_k: 95.0,
                risk_free_rate_r: r,
                effective_volatility_sigma_bar: vol,
                fast_skew_parameter_v3_eps: -0.0012 - credit_bump * 0.05,
                fast_convexity_parameter_v2_eps: 0.0018,
                slow_skew_parameter_v1_delta: -0.0015,
                slow_level_parameter_v0_delta: 0.0010,
                maturity_years_t: 0.75,
            },
        );
        let ga_ftd = ql_xva::evaluate_gregory_albanese_bilateral_cva_dva_first_to_default_survival(
            &ql_xva::GregoryAlbaneseBilateralFtdCvaDvaSpec {
                expected_positive_exposure_epe: 15.0,
                expected_negative_exposure_ene: 10.0,
                counterparty_hazard_rate_lambda_c: 0.030 + credit_bump,
                own_bank_hazard_rate_lambda_b: 0.018 + credit_bump * 0.5,
                counterparty_lgd_c: 0.60,
                own_bank_lgd_b: 0.60,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let ew_vpp = ql_commodities::evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip(
            &ql_commodities::EydelandWolyniecVppHeatRateStripSpec {
                power_forward_price_f_p: 68.0 + inflation * 50.0,
                gas_forward_price_f_g: 7.50,
                plant_heat_rate_h: 7.80,
                variable_om_cost_v_om: 3.50,
                startup_cost_per_mwh_c_start: 2.20,
                power_volatility_sigma_p: 0.48 + credit_bump,
                gas_volatility_sigma_g: 0.36,
                power_gas_correlation_rho: 0.72,
                risk_free_rate_r: r,
                maturity_years_t: 0.75,
            },
        );
        let fx_corr = ql_fx::evaluate_wystup_fx_corridor_accrual_range_note(
            &ql_fx::WystupFxCorridorAccrualRangeNoteSpec {
                notional_coupon_max: 100.0,
                spot_fx_s0: 1.10,
                lower_corridor_barrier_l: 1.04,
                upper_corridor_barrier_u: 1.18,
                domestic_rate_rd: r + 0.005,
                foreign_rate_rf: 0.030,
                atm_volatility_sigma: 0.095 + credit_bump,
                log_moneyness_skew_slope_xi: -0.12,
                maturity_years_t: 1.0,
            },
        );
        let rr_do = ql_exotics::evaluate_rubinstein_reiner_down_and_out_call_with_rebate(
            &ql_exotics::RubinsteinReinerDownAndOutCallRebateSpec {
                spot_price_s0: spot,
                strike_price_k: 100.0,
                lower_barrier_h: 90.0,
                cash_rebate_at_hit_r_rebate: 3.0,
                risk_free_rate_r: r,
                dividend_yield_q: 0.02,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let lns_gls = ql_portfolio::evaluate_lewellen_nagel_shanken_gls_cross_sectional_r2(
            &ql_portfolio::LewellenNagelShankenGlsCrossSectionalSpec {
                factor_tangency_sharpe_sq_theta_f2: 0.085,
                residual_alpha_quadratic_form_q_alpha: 0.140 + credit_bump,
                cross_sectional_ols_r2: 0.78,
                factor_risk_premium_lambda: 0.0065,
                factor_variance_sigma_f2: 0.0022,
                sample_length_t: 240.0,
                test_asset_count_n: 25.0,
                factor_count_k: 3.0,
            },
        );
        let bc_exp = ql_portfolio::evaluate_bellini_cont_expectile_evar_influence_sensitivity(
            &ql_portfolio::BelliniContExpectileRiskMeasureSpec {
                loss_mean_mu: 0.50 + credit_bump * 10.0,
                loss_volatility_sigma: 2.20,
                tail_asymmetry_tau: 0.975,
                contamination_loss_shock_z: 8.50,
                contamination_weight_eps: 0.02,
            },
        );
        let bns_rq = ql_dml::evaluate_barndorff_nielsen_shephard_realized_quarticity_iv_ci(
            &ql_dml::BarndorffNielsenShephardRealizedQuarticitySpec {
                realized_variance_rv: 0.032 + credit_bump,
                realized_quarticity_rq: 0.0014,
                intraday_sampling_frequency_m: 78.0,
                first_order_autocovariance_gamma1: -0.0008,
                normal_critical_z_alpha: 1.96,
            },
        );
        let km_cyc = ql_models::evaluate_kiyotaki_moore_collateral_credit_cycle_multiplier(
            &ql_models::KiyotakiMooreCollateralCycleSpec {
                land_asset_price_q: 100.0 - credit_bump * 100.0,
                expected_next_land_price_q_next: 101.0,
                gross_interest_rate_r_gross: 1.0 + r,
                farmer_productivity_a: 6.0,
                existing_land_holding_k_prev: 10.0,
                inherited_debt_b_prev: 880.0,
                intertemporal_price_elasticity_eta: 0.02,
            },
        );

        lm_cap.backward_looking_in_arrears_caplet_npv
            + cpi_seas.zciis_seasonal_mark_to_market_npv
            + ba_cir.correlated_ssrd_cds_fair_spread_bps * 0.01
            + fpss_sv.multiscale_sv_call_option_price * 0.1
            + ga_ftd.net_bilateral_bcva
            + ew_vpp.vpp_startup_adjusted_tolling_value * 0.1
            + fx_corr.smile_adjusted_range_note_pv * 0.1
            + rr_do.total_down_and_out_call_with_rebate * 0.1
            + lns_gls.gls_cross_sectional_r2 * 10.0
            + bc_exp.contaminated_robust_expectile
            + bns_rq.log_transformed_iv_upper_bound * 100.0
            + km_cyc.dynamic_credit_cycle_amplification_multiplier
    }
    pub fn lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let jsz_gdtsm = ql_curves::evaluate_joslin_singleton_zhu_canonical_gdtsm_term_premium(
            &ql_curves::JoslinSingletonZhuGdtsmSpec {
                risk_neutral_long_run_mean_r_inf_q: r + 0.015,
                physical_long_run_mean_r_inf_p: 0.032 + inflation * 0.1,
                risk_neutral_decay_speeds_kappa_q: [0.06, 0.45, 1.20],
                physical_decay_speeds_kappa_p: [0.18, 0.65, 1.45],
                canonical_factor_states_x: [0.015 + credit_bump, 0.008, 0.012],
                factor_volatilities_sigma: [0.0065, 0.0080, 0.0095],
                bond_maturity_years_tau: 10.0,
            },
        );
        let tf_cb = ql_fixed_income::evaluate_tsiveriotis_fernandes_convertible_bond_credit_split(
            &ql_fixed_income::TsiveriotisFernandesConvertibleSpec {
                spot_price_s0: spot,
                conversion_ratio_gamma: 1.0,
                par_value: 100.0,
                coupon_rate_c: 0.035,
                risk_free_rate_r: r,
                credit_spread_sc: 0.025 + credit_bump,
                dividend_yield_q: 0.015,
                equity_volatility_sigma: vol,
                call_price_k: 112.0,
                maturity_years_t: 2.0,
            },
        );
        let lg_cdo = ql_credit::evaluate_laurent_gregory_fft_recursion_synthetic_cdo_tranche(
            &ql_credit::LaurentGregoryRecursiveCdoSpec {
                attachment_point_ka: 0.03,
                detachment_point_kd: 0.07,
                recovery_rate_rec: 0.40,
                tier_hazard_rates_lambda: [
                    0.008 + credit_bump,
                    0.020 + credit_bump,
                    0.045 + credit_bump,
                ],
                tier_factor_correlations_rho: [0.25, 0.35, 0.45],
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let bj_vs = ql_vol::evaluate_broadie_jain_discrete_sampling_jump_variance_swap(
            &ql_vol::BroadieJainVarianceSwapSpec {
                diffusive_volatility_sigma: vol,
                jump_intensity_lambda_j: 0.60 + credit_bump * 5.0,
                mean_log_jump_mu_j: -0.08,
                std_log_jump_delta_j: 0.12,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                discrete_sampling_periods_n: 52.0,
                maturity_years_t: 1.0,
            },
        );
        let ac_kva = ql_xva::evaluate_albanese_crepey_sa_ccr_kadd_capital_xva(
            &ql_xva::AlbaneseCrepeySaCcrKvaSpec {
                trade_notional: 1_000.0,
                current_mtm_value_v0: 8.0 + credit_bump * 50.0,
                collateral_balance_c0: 12.0,
                counterparty_risk_weight_w: 0.50,
                capital_hurdle_rate_h: 0.10,
                ois_discount_rate_r: r,
                joint_hazard_rate_lambda: 0.025 + credit_bump,
                maturity_years_t: 5.0,
            },
        );
        let rss_stor =
            ql_commodities::evaluate_routledge_seppi_spatt_equilibrium_storage_forward_curve(
                &ql_commodities::RoutledgeSeppiSpattStorageSpec {
                    spot_price_s0: 78.0 + inflation * 50.0,
                    physical_inventory_buffer_i0: 0.12,
                    net_demand_shock_x0: 0.18 + credit_bump,
                    demand_mean_reversion_kappa_x: 1.40,
                    demand_shock_volatility_sigma_x: 0.26,
                    storage_cost_rate_cs: 0.035,
                    risk_free_rate_r: r,
                    forward_maturity_years_tau: 0.75,
                },
            );
        let dtv_cip = ql_fx::evaluate_du_tepper_verdelhan_cip_treasury_basis_arb(
            &ql_fx::DuTepperVerdelhanCipBasisSpec {
                spot_fx_s0: 1.0800,
                outright_forward_fx_f0: 1.1015,
                domestic_usd_ois_rate_rd: r + 0.01,
                foreign_ois_rate_rf: 0.0310,
                domestic_ust_yield: r + 0.0085,
                foreign_gov_yield: 0.0325,
                dealer_slr_balance_sheet_cost: 0.0018 + credit_bump * 0.2,
                dealer_funding_spread_delta_s: 0.0012,
                maturity_years_tau: 1.0,
            },
        );
        let wc_cliq = ql_exotics::evaluate_wilmott_carr_globally_floored_locally_capped_cliquet(
            &ql_exotics::WilmottCarrCliquetSpec {
                notional_n: 100.0,
                local_return_floor_f_loc: 0.00,
                local_return_cap_c_loc: 0.08,
                global_cumulative_floor_f_glob: 0.10,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                forward_atm_volatility_sigma: vol,
                forward_skew_slope_eta: -0.25 - credit_bump,
                maturity_years_t: 2.0,
            },
        );
        let bsv_ppp = ql_portfolio::evaluate_brandt_santa_clara_valkanov_crra_characteristic_tilt(
            &ql_portfolio::BrandtSantaClaraValkanovCharacteristicTiltSpec {
                benchmark_expected_return_mu_b: 0.068 - credit_bump,
                benchmark_volatility_sigma_b: vol,
                value_characteristic_premium_mu_v: 0.042,
                momentum_characteristic_premium_mu_m: 0.058,
                value_characteristic_vol_sigma_v: 0.115 + credit_bump,
                momentum_characteristic_vol_sigma_m: 0.145,
                characteristic_correlation_rho_vm: -0.35,
                relative_risk_aversion_gamma: 5.0,
            },
        );
        let aj_evar = ql_portfolio::evaluate_ahmadi_javid_follmer_schied_entropic_var_evar(
            &ql_portfolio::AhmadiJavidEntropicVarSpec {
                mean_diffusive_loss_mu: 0.25 + credit_bump * 5.0,
                diffusive_loss_volatility_sigma: 1.50,
                jump_loss_intensity_lambda: 0.15 + credit_bump,
                mean_jump_loss_mu_j: 2.20,
                jump_loss_volatility_delta_j: 0.85,
                tail_significance_alpha: 0.025,
            },
        );
        let ek_bekk = ql_dml::evaluate_engle_kroner_bekk_bivariate_garch_covariance(
            &ql_dml::EngleKronerBekkGarchSpec {
                cholesky_intercept_c: [0.04 + credit_bump, 0.015, 0.035],
                arch_spillover_matrix_a: [0.28, 0.08, 0.06, 0.25],
                garch_persistence_diag_b: [0.92, 0.90],
                lagged_innovations_eps: [-0.025, -0.018],
                lagged_covariance_h: [0.032, 0.018, 0.028],
            },
        );
        let gab_dis = ql_models::evaluate_gabaix_variable_rare_disaster_resilience_equity_puzzle(
            &ql_models::GabaixVariableRareDisasterSpec {
                risk_aversion_gamma: 4.0,
                time_preference_rho: 0.035,
                normal_consumption_growth_gc: 0.020,
                consumption_volatility_sigma_c: 0.018,
                current_disaster_probability_pt: 0.022 + credit_bump * 0.5,
                baseline_disaster_probability_p_star: 0.017,
                consumption_disaster_recovery_bc: 0.72,
                dividend_disaster_recovery_fd: 0.55,
                resilience_mean_reversion_phi_h: 0.18,
            },
        );

        jsz_gdtsm.canonical_yield_term_premium * 100.0
            + tf_cb.convertible_bond_value * 0.1
            + lg_cdo.fair_tranche_running_spread_bps * 0.01
            + bj_vs.convexity_adjusted_vol_swap_strike * 10.0
            + ac_kva.lifetime_sa_ccr_kva
            + rss_stor.equilibrium_forward_price * 0.1
            + dtv_cip.net_actionable_cip_arbitrage_alpha_bps * 0.1
            + wc_cliq.total_cliquet_present_value * 0.1
            + bsv_ppp.crra_certainty_equivalent_return * 100.0
            + aj_evar.jump_diffusion_entropic_var
            + ek_bekk.conditional_correlation_rho12 * 10.0
            + gab_dis.equilibrium_price_dividend_ratio
    }
    pub fn jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let bian_q = ql_curves::evaluate_bianchetti_two_curve_forward_measure_quanto_convexity(
            &ql_curves::BianchettiTwoCurveForwardQuantoSpec {
                notional_n: 1_000.0,
                projection_forward_rate_f0: r + 0.012,
                tenor_basis_spread_s0: 0.0035 + credit_bump * 0.5,
                caplet_strike_k: 0.042,
                fixing_expiry_years_t1: 2.0,
                accrual_fraction_tau: 0.50,
                forward_rate_vol_sigma_f: 0.25 + credit_bump,
                basis_spread_vol_sigma_s: 0.35,
                forward_basis_correlation_rho: 0.65,
                ois_discount_rate_r: r,
            },
        );
        let bdt_cb = ql_fixed_income::evaluate_black_derman_toy_lognormal_short_rate_callable_bond(
            &ql_fixed_income::BlackDermanToyCallableBondSpec {
                par_value: 100.0,
                coupon_rate_c: 0.055,
                initial_short_rate_r0: r,
                term_slope_s: 0.002,
                short_rate_log_vol_sigma: 0.22 + credit_bump,
                option_adjusted_spread_oas: 0.006 + credit_bump,
                call_price_k: 101.0,
                maturity_years_t: 4.0,
            },
        );
        let sj_swp = ql_credit::evaluate_schonbucher_jamshidian_survival_measure_payer_cds_swaption(
            &ql_credit::SchonbucherJamshidianCdsSwaptionSpec {
                notional_n: 1_000.0,
                short_term_hazard_rate_lambda0: 0.018 + credit_bump,
                forward_hazard_rate_lambda_fwd: 0.022 + credit_bump,
                strike_spread_k: 0.012,
                spread_volatility_sigma_s: 0.45,
                recovery_rate_rec: 0.40,
                risk_free_rate_r: r,
                option_expiry_years_ta: 1.0,
                cds_maturity_years_tb: 6.0,
            },
        );
        let bates_svj = ql_vol::evaluate_bates_svj_stochastic_vol_jump_characteristic_smile(
            &ql_vol::BatesSvjSmileSpec {
                spot_price_s0: spot,
                strike_price_k: 93.0,
                risk_free_rate_r: r,
                initial_variance_v0: vol * vol,
                long_run_variance_theta: 0.045,
                mean_reversion_kappa: 1.80,
                vol_of_vol_xi: 0.38,
                spot_vol_correlation_rho: -0.65,
                jump_intensity_lambda_j: 0.35 + credit_bump * 5.0,
                mean_log_jump_mu_j: -0.07,
                std_log_jump_delta_j: 0.11,
                maturity_years_t: 0.75,
            },
        );
        let bp_ccp = ql_xva::evaluate_brigo_pallavicini_ccp_default_fund_cva_mva(
            &ql_xva::BrigoPallaviciniCcpXvaSpec {
                portfolio_mtm_volatility_sigma_v: 45.0,
                margin_period_of_risk_years: 10.0 / 252.0,
                im_confidence_alpha: 0.99,
                im_funding_spread_s_im: 0.0085 + credit_bump * 0.5,
                clearing_member_hazard_lambda_cm: 0.015 + credit_bump,
                default_fund_scaling_theta_df: 0.18,
                default_fund_mutualization_spread_s_df: 0.025,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let dmt_crk =
            ql_commodities::evaluate_dempster_medova_tang_cointegrated_crack_spread_option(
                &ql_commodities::DempsterMedovaTangCrackSpreadSpec {
                    refined_product_forward_f1: 96.0 + inflation * 50.0,
                    crude_feedstock_forward_f2: 78.0,
                    refining_strike_k: 14.0,
                    product_volatility_sigma1: 0.34 + credit_bump,
                    crude_volatility_sigma2: 0.30,
                    product_crude_correlation_rho: 0.82,
                    cointegration_reversion_speed_kappa: 1.50,
                    risk_free_rate_r: r,
                    maturity_years_t: 1.5,
                },
            );
        let pit_prdc =
            ql_fx::evaluate_piterbarg_fx_power_reverse_dual_currency_prdc_bermudan_redemption(
                &ql_fx::PiterbargPrdcBermudanSpec {
                    notional_n: 1_000.0,
                    spot_fx_s0: 105.0,
                    foreign_coupon_gearing_gf: 0.10,
                    domestic_funding_deduction_gd: 0.06,
                    domestic_rate_rd: r * 0.5,
                    foreign_rate_rf: 0.045,
                    fx_volatility_sigma_fx: 0.14 + credit_bump,
                    domestic_rate_vol_sigma_d: 0.012,
                    domestic_fx_correlation_rho: 0.45,
                    target_redemption_cap_c_target: 0.12,
                    maturity_years_t: 4.0,
                },
            );
        let bgk_bar =
            ql_exotics::evaluate_broadie_glasserman_kou_discrete_barrier_continuity_correction(
                &ql_exotics::BroadieGlassermanKouDiscreteBarrierSpec {
                    spot_price_s0: spot,
                    strike_price_k: 100.0,
                    contractual_down_barrier_h: 91.0,
                    discrete_monitoring_dates_m: 52.0,
                    risk_free_rate_r: r,
                    dividend_yield_q: 0.02,
                    volatility_sigma: vol,
                    maturity_years_t: 1.0,
                },
            );
        let ck_apt = ql_portfolio::evaluate_connor_korajczyk_asymptotic_pca_apt_pricing_error(
            &ql_portfolio::ConnorKorajczykAptSpec {
                asset_expected_excess_return_mu: 0.088,
                asset_pca_factor_betas: [1.15, 0.45],
                pca_factor_risk_premia_lambda: [0.055 - credit_bump, 0.032],
                pca_factor_variances: [0.025, 0.014],
                idiosyncratic_residual_variance_sigma_e2: 0.036,
                cross_sectional_asset_count_n: 250.0,
                time_series_sample_length_t: 180.0,
            },
        );
        let afs_exec =
            ql_portfolio::evaluate_alfonsi_fruth_schied_bucket_transient_impact_execution(
                &ql_portfolio::AlfonsiFruthSchiedExecutionSpec {
                    parent_order_size_x0: 100.0,
                    discrete_intervals_n: 10.0,
                    lob_resilience_speed_rho: 2.50,
                    transient_impact_coeff_eta: 0.04 + credit_bump,
                    permanent_impact_coeff_gamma: 0.01,
                    asset_volatility_sigma: vol,
                    inventory_risk_aversion_lambda: 0.05,
                    horizon_t: 1.0,
                },
            );
        let ham_ms = ql_dml::evaluate_hamilton_markov_switching_regime_filter_step(
            &ql_dml::HamiltonMarkovSwitchingSpec {
                observed_realization_yt: -0.018 - credit_bump,
                expansion_mean_mu1: 0.008,
                recession_mean_mu2: -0.022,
                expansion_volatility_sigma1: 0.012,
                recession_volatility_sigma2: 0.028,
                stay_expansion_prob_p11: 0.95,
                stay_recession_prob_p22: 0.85,
                prior_expansion_prob_xi1: 0.80,
            },
        );
        let by_lrr = ql_models::evaluate_bansal_yaron_long_run_risk_stochastic_vol_equity_premium(
            &ql_models::BansalYaronLongRunRiskSpec {
                risk_aversion_gamma: 7.5,
                eis_psi: 1.5,
                growth_persistence_rho_x: 0.975,
                growth_shock_scale_phi_e: 0.044,
                volatility_persistence_nu: 0.985,
                vol_of_vol_sigma_w: 0.00025,
                current_consumption_vol_sigma_t: 0.027 + credit_bump * 0.5,
                dividend_leverage_phi_d: 3.0,
                campbell_shiller_kappa1: 0.99,
            },
        );

        bian_q.bianchetti_quanto_adjusted_caplet_pv
            + bdt_cb.bdt_callable_bond_price * 0.1
            + sj_swp.full_non_knockout_swaption_pv * 0.1
            + bates_svj.bates_svj_call_option_price * 0.1
            + bp_ccp.total_ccp_clearing_xva
            + dmt_crk.cointegrated_crack_spread_option_price * 0.1
            + pit_prdc.tarn_capped_prdc_present_value * 0.1
            + bgk_bar.bgk_discrete_down_and_out_call_price * 0.1
            + ck_apt.connor_korajczyk_wald_statistic
            + afs_exec.optimal_afs_execution_cost * 0.1
            + ham_ms.next_step_conditional_forecast * 100.0
            + by_lrr.conditional_equity_risk_premium * 100.0
    }
    pub fn bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let ap_cms = ql_curves::evaluate_andersen_piterbarg_lsrm_cms_spread_option(
            &ql_curves::AndersenPiterbargCmsSpreadSpec {
                cms1_forward_swap_rate_s1: r + 0.010,
                cms1_tenor_years_tau1: 10.0,
                cms1_normal_volatility_sigma1: 0.0085 + credit_bump * 0.1,
                cms2_forward_swap_rate_s2: r + 0.002,
                cms2_tenor_years_tau2: 2.0,
                cms2_normal_volatility_sigma2: 0.0095,
                swap_rate_correlation_rho12: 0.78,
                leg1_gearing_w1: 1.0,
                leg2_gearing_w2: 1.0,
                spread_strike_k: 0.005,
                mean_reversion_kappa: 0.08,
                ois_discount_rate_r: r,
                fixing_horizon_years_t0: 3.0,
            },
        );
        let fab_z = ql_fixed_income::evaluate_fabozzi_cmo_accrual_z_bond_sequential_lockout_tranche(
            &ql_fixed_income::FabozziCmoAccrualZBondSpec {
                senior_tranche_initial_balance_ba: 75.0,
                z_bond_initial_balance_bz: 25.0,
                tranche_coupon_rate_c: 0.055,
                conditional_prepayment_rate_cpr: 0.16 - credit_bump,
                scheduled_amortization_rate: 0.12,
                base_discount_yield_y: r,
                z_bond_oas_spread: 0.012 + credit_bump,
            },
        );
        let dg_cdo = ql_credit::evaluate_duffie_garleanu_affine_jump_intensity_cdo_tranche(
            &ql_credit::DuffieGarleanuAffineJumpCdoSpec {
                initial_common_intensity_lambda_c0: 0.015 + credit_bump,
                common_reversion_speed_kappa_c: 0.45,
                common_long_run_theta_c: 0.018,
                common_feller_volatility_sigma_c: 0.09,
                common_jump_arrival_rate_mu_j: 0.25 + credit_bump * 2.0,
                common_mean_jump_size_eta_j: 0.03,
                idiosyncratic_intensity_lambda_id: 0.012 + credit_bump,
                loss_given_default_lgd: 0.60,
                attachment_point_ka: 0.03,
                detachment_point_kd: 0.07,
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let cl_vol = ql_vol::evaluate_carr_lee_robust_volatility_swap_convexity_correction(
            &ql_vol::CarrLeeVolatilitySwapSpec {
                initial_variance_v0: vol * vol,
                long_run_variance_theta: 0.040,
                mean_reversion_speed_kappa: 1.80,
                vol_of_vol_xi: 0.45 + credit_bump,
                realized_variance_skewness_gamma: 0.65,
                maturity_years_t: 1.0,
            },
        );
        let gr_kva = ql_xva::evaluate_green_basel_leverage_ratio_rwa_hurdle_kva(
            &ql_xva::GreenBaselLeverageRatioKvaSpec {
                initial_sa_ccr_ead0: 18.0 + credit_bump * 50.0,
                gross_trade_notional: 500.0,
                sa_ccr_addon_factor_psi: 0.04,
                counterparty_risk_weight_w: 0.50,
                cet1_capital_ratio_k: 0.105,
                leverage_ratio_requirement_gamma: 0.035,
                cva_capital_loading_beta: 0.08,
                capital_hurdle_rate_h: 0.11 + credit_bump,
                capital_funding_offset_phi: 0.35,
                counterparty_hazard_rate_lambda: 0.015 + credit_bump,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let cs_oil = ql_commodities::evaluate_cortazar_schwartz_three_factor_oil_futures_curve(
            &ql_commodities::CortazarSchwartzThreeFactorSpec {
                spot_price_s0: 78.0 + inflation * 50.0,
                initial_convenience_yield_delta0: 0.08,
                initial_long_term_drift_v0: 0.035,
                risk_neutral_convenience_mean_alpha2: 0.05,
                risk_neutral_long_run_drift_vbar: 0.03,
                convenience_reversion_speed_kappa2: 1.40,
                drift_reversion_speed_kappa3: 0.35,
                spot_volatility_sigma1: 0.32 + credit_bump,
                convenience_volatility_sigma2: 0.24,
                drift_volatility_sigma3: 0.08,
                spot_convenience_correlation_rho12: 0.72,
                maturity_years_tau: 2.5,
            },
        );
        let fst_xccs = ql_fx::evaluate_fujii_shimada_takahashi_collateralized_xccs_curve_discount(
            &ql_fx::FujiiShimadaTakahashiXccsSpec {
                spot_fx_s0: 1.08,
                domestic_ois_rate_rd: r,
                foreign_ois_rate_rf: 0.030,
                xccy_collateral_basis_spread_b: -0.0025 - credit_bump * 0.5,
                foreign_swap_coupon_cf: 0.032,
                foreign_rate_volatility_sigma_f: 0.009,
                fx_volatility_sigma_fx: 0.11 + credit_bump,
                foreign_fx_correlation_rho: 0.35,
                maturity_years_t: 4.0,
            },
        );
        let ph_tsw = ql_exotics::evaluate_pechtl_haug_discrete_time_switch_accumulator_option(
            &ql_exotics::PechtlHaugTimeSwitchSpec {
                spot_price_s0: spot,
                strike_price_k: 100.0,
                accrual_rate_a: 25.0,
                accumulated_occupation_time_t_acc: 0.25,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                smile_skew_slope_nu: -0.08 - credit_bump,
                maturity_years_t: 1.5,
            },
        );
        let bcf_mvar =
            ql_portfolio::evaluate_boudt_cornish_fisher_modified_var_risk_budget_allocation(
                &ql_portfolio::BoudtCornishFisherMvarSpec {
                    initial_weight_w1: 0.60,
                    target_risk_budget_b1: 0.50,
                    volatility_sigma1: vol,
                    volatility_sigma2: 0.14,
                    return_correlation_rho: 0.35,
                    loss_skewness_loading_s1: 0.45 + credit_bump * 2.0,
                    loss_skewness_loading_s2: 0.15,
                    excess_kurtosis_loading_k1: 1.80,
                    excess_kurtosis_loading_k2: 0.60,
                    normal_tail_quantile_z: 2.3263,
                },
            );
        let ru_quad = ql_portfolio::evaluate_rockafellar_uryasev_mixed_cvar_risk_quadrangle(
            &ql_portfolio::RockafellarUryasevMixedCvarSpec {
                mean_loss_mu: 0.012 + credit_bump,
                loss_volatility_sigma: 0.045 + credit_bump,
                moderate_confidence_alpha1: 0.95,
                extreme_confidence_alpha2: 0.99,
                moderate_spectral_weight_lambda1: 0.60,
                tail_kurtosis_amplification_kappa: 1.50,
            },
        );
        let hgg_is = ql_dml::evaluate_hasbrouck_gonzalo_granger_price_discovery_information_share(
            &ql_dml::HasbrouckGonzaloGrangerSpec {
                venue1_error_correction_speed_alpha1: 0.08 + credit_bump,
                venue2_error_correction_speed_alpha2: 0.24,
                venue1_innovation_volatility_sigma1: 0.015,
                venue2_innovation_volatility_sigma2: 0.018,
                cross_venue_innovation_correlation_rho: 0.55,
            },
        );
        let hk_sdf = ql_models::evaluate_he_krishnamurthy_intermediary_capital_constraint_sdf(
            &ql_models::HeKrishnamurthyIntermediarySdfSpec {
                intermediary_wealth_share_wi: (0.14 - credit_bump).max(0.05),
                capital_constraint_threshold_w_star: 0.28,
                equity_intermediation_pledgeability_m: 0.35,
                specialist_risk_aversion_gamma: 2.20,
                fundamental_growth_rate_ga: 0.020,
                asset_return_volatility_sigma_a: vol,
                subjective_discount_rate_rho: 0.025,
            },
        );

        ap_cms.cms_spread_call_present_value * 100.0
            + fab_z.z_bond_present_value * 0.1
            + dg_cdo.tranche_expected_loss_pv * 10.0
            + cl_vol.volatility_swap_fair_strike_k_vol * 10.0
            + gr_kva.total_dual_constraint_kva
            + cs_oil.three_factor_futures_price * 0.1
            + fst_xccs.collateralized_xccs_present_value * 10.0
            + ph_tsw.time_switch_option_present_value
            + bcf_mvar.tail_risk_parity_weight_w1 * 10.0
            + ru_quad.quadrangle_mixed_cvar_risk * 10.0
            + hgg_is.hasbrouck_midpoint_information_share * 10.0
            + hk_sdf.conditional_intermediary_equity_premium * 100.0
    }
    pub fn andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let pit_ctd = ql_curves::evaluate_piterbarg_multi_currency_ctd_collateral_choice_discount(
            &ql_curves::PiterbargMultiCurrencyCtdCollateralSpec {
                domestic_ois_rate_rd: r,
                domestic_collateral_spread_s1: 0.0010,
                foreign_synthetic_collateral_spread_s2: 0.0025 + credit_bump * 0.5,
                collateral_spread_diff_volatility_sigma_x: 0.0035 + credit_bump * 0.1,
                spread_reversion_speed_kappa_x: 0.45,
                maturity_years_t: 5.0,
            },
        );
        let tuck_mw = ql_fixed_income::evaluate_tuckman_corporate_make_whole_and_par_call_bond(
            &ql_fixed_income::TuckmanMakeWholeCallableBondSpec {
                par_value_f: 100.0,
                coupon_rate_c: 0.056,
                treasury_benchmark_yield_y: r,
                issuer_credit_spread_s_corp: 0.014 + credit_bump,
                make_whole_treasury_spread_s_mw: 0.0030,
                credit_spread_volatility_sigma_s: 0.0065 + credit_bump * 0.1,
                treasury_yield_volatility_sigma_y: 0.0085,
                make_whole_horizon_t_mw: 2.0,
                par_call_window_start_t_par: 4.5,
                maturity_years_t_mat: 5.0,
            },
        );
        let bm_cmcds = ql_credit::evaluate_brigo_morini_constant_maturity_cds_convexity(
            &ql_credit::BrigoMoriniCmcdsSpec {
                notional_n: 1_000.0,
                forward_reference_cds_spread_s0: 0.015 + credit_bump,
                reference_cds_tenor_tau_cm: 5.0,
                cds_spread_volatility_sigma_s: 0.45,
                survival_hazard_rate_lambda: 0.025 + credit_bump,
                loss_given_default_lgd: 0.60,
                cmcds_cap_strike_k: 0.025,
                risk_free_rate_r: r,
                contract_maturity_years_t: 4.0,
            },
        );
        let hag_sabr = ql_vol::evaluate_hagan_shifted_sabr_negative_rate_implied_volatility(
            &ql_vol::HaganShiftedSabrSpec {
                forward_rate_f: -0.004 + credit_bump * 0.5,
                strike_rate_k: -0.008,
                displacement_shift_s: 0.030,
                sabr_alpha: 0.045,
                sabr_beta: 0.50,
                sabr_vol_of_vol_nu: 0.40 + credit_bump,
                sabr_correlation_rho: -0.35,
                discount_rate_r: r,
                expiry_years_t: 2.0,
            },
        );
        let gk_tva = ql_xva::evaluate_green_kenyon_tax_valuation_adjustment_tva(
            &ql_xva::GreenKenyonTvaSpec {
                initial_epe0: 25.0 + credit_bump * 50.0,
                initial_ene0: 14.0,
                counterparty_hazard_rate_lambda_c: 0.022 + credit_bump,
                bank_hazard_rate_lambda_b: 0.012,
                loss_given_default_lgd: 0.60,
                unsecured_funding_spread_sf: 0.0095 + credit_bump,
                corporate_tax_rate_tau: 0.25,
                funding_interest_deductibility_phi: 0.85,
                dva_tax_realizability_psi: 0.40,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let ew_inter =
            ql_commodities::evaluate_eydeland_wolyniec_interconnector_transmission_spread_option(
                &ql_commodities::EydelandWolyniecInterconnectorSpec {
                    sink_pool_forward_f_sink: 88.0 + inflation * 50.0,
                    source_pool_forward_f_source: 68.0,
                    transmission_line_loss_fraction: 0.045,
                    wheeling_tariff_strike_k: 6.0,
                    sink_volatility_sigma1: 0.42 + credit_bump,
                    source_volatility_sigma2: 0.36,
                    cross_pool_correlation_rho: 0.74,
                    forced_outage_hazard_lambda: 0.80,
                    repair_rate_mu: 8.50,
                    risk_free_rate_r: r,
                    maturity_years_t: 1.0,
                },
            );
        let wc_str = ql_fx::evaluate_wystup_clark_fx_smile_strangle_market_convention(
            &ql_fx::WystupClarkFxSmileStrangleSpec {
                spot_fx_s0: 1.10,
                domestic_rate_rd: r,
                foreign_rate_rf: 0.030,
                atm_volatility_sigma_atm: 0.095 + credit_bump,
                risk_reversal_25d_rr: -0.015,
                butterfly_25d_bf: 0.0035,
                expiry_years_t: 1.0,
            },
        );
        let ki_bar = ql_exotics::evaluate_kunitomo_ikeda_curved_exponential_barrier_option(
            &ql_exotics::KunitomoIkedaCurvedBarrierSpec {
                spot_price_s0: spot,
                strike_price_k: 100.0,
                initial_down_barrier_h0: 88.0,
                barrier_curvature_rate_delta_h: 0.02,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let clmx_var =
            ql_portfolio::evaluate_campbell_lettau_malkiel_xu_hierarchical_variance_decomposition(
                &ql_portfolio::CampbellLettauMalkielXuSpec {
                    market_volatility_sigma_m: vol,
                    industry_residual_volatility_sigma_ind: 0.12 + credit_bump,
                    firm_idiosyncratic_volatility_sigma_firm: 0.32,
                    industry_herfindahl_index_h_ind: 0.08,
                    portfolio_stock_count_n: 30.0,
                },
            );
        let kyle_eq = ql_portfolio::evaluate_kyle_multi_period_insider_market_depth_equilibrium(
            &ql_portfolio::KyleMultiPeriodEquilibriumSpec {
                prior_fundamental_variance_sigma0: 4.0 + credit_bump * 10.0,
                noise_trader_volatility_sigma_u: 12.0,
                insider_valuation_wedge_v_minus_p0: 2.5,
                auction_step_dt: 0.5,
            },
        );
        let hl_rv = ql_dml::evaluate_hansen_lunde_overnight_intraday_realized_variance_combination(
            &ql_dml::HansenLundeOvernightRvSpec {
                realized_overnight_sq_return_rv1: 0.00012 + credit_bump * 0.001,
                realized_intraday_variance_rv2: 0.00038,
                unconditional_overnight_mean_mu1: 0.00010,
                unconditional_intraday_mean_mu2: 0.00035,
                unconditional_full_day_iv_mu0: 0.00042,
                overnight_estimator_variance_eta1_sq: 2.5e-8,
                intraday_estimator_variance_eta2_sq: 1.8e-8,
                estimator_covariance_eta12: 0.4e-8,
            },
        );
        let dd_gp = ql_models::evaluate_diamond_dybvig_goldstein_pauzner_bank_run_global_game(
            &ql_models::DiamondDybvigGoldsteinPauznerBankRunSpec {
                short_term_deposit_payout_r1: 1.06,
                long_term_asset_return_r_bar: 1.55,
                liquid_reserve_ratio_l: (0.22 - credit_bump).max(0.05),
                fire_sale_recovery_fraction_psi: 0.65,
                impatient_depositor_fraction_lambda: 0.30,
                fundamental_prior_mean_mu_theta: 0.85 - credit_bump,
                fundamental_volatility_sigma_theta: 0.18,
            },
        );

        pit_ctd.ctd_collateral_discount_factor * 10.0
            + tuck_mw.callable_corporate_bond_price * 0.1
            + bm_cmcds.capped_cmcds_floating_leg_pv * 0.1
            + hag_sabr.shifted_bachelier_call_price * 100.0
            + gk_tva.post_tax_net_xva_with_tva
            + ew_inter.expected_interconnector_congestion_rent_pv
            + wc_str.smile_adjusted_strangle_present_value * 10.0
            + ki_bar.kunitomo_ikeda_down_and_out_call_price * 0.1
            + clmx_var.finite_n_portfolio_volatility * 10.0
            + kyle_eq.expected_insider_information_profit
            + hl_rv.hansen_lunde_combined_realized_volatility * 100.0
            + dd_gp.expected_depositor_welfare_value * 10.0
    }
    pub fn piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let merc_fva = ql_curves::evaluate_mercurio_post_crisis_fva_adjusted_multi_curve_swap_rate(
            &ql_curves::MercurioFvaMultiCurveSwapSpec {
                ois_short_rate_r0: r,
                ois_curve_slope_theta: 0.004,
                initial_multiplicative_ibor_basis_b0: 0.0025 + credit_bump * 0.5,
                basis_term_decay_kappa_b: 0.25,
                funding_borrowing_spread_s_f: 0.0065 + credit_bump,
                swap_rate_volatility_sigma_s: 0.0085,
                swap_maturity_years_t: 4.0,
            },
        );
        let ing_cb =
            ql_fixed_income::evaluate_ingersoll_contingent_call_forced_conversion_convertible_bond(
                &ql_fixed_income::IngersollConvertibleBondSpec {
                    stock_price_s0: spot,
                    conversion_ratio_cr: 1.0,
                    straight_bond_coupon_c: 0.04,
                    call_trigger_price_k_call: 125.0,
                    call_protection_years_t_call: 2.0,
                    maturity_years_t: 5.0,
                    credit_adjusted_discount_rate_r_c: r + 0.02 + credit_bump,
                    equity_volatility_sigma: vol,
                },
            );
        let hw_cdo = ql_credit::evaluate_hull_white_implied_copula_bespoke_cdo_tranche_calibration(
            &ql_credit::HullWhiteImpliedCopulaCdoSpec {
                base_hazard_rate_lambda0: 0.015 + credit_bump,
                recession_hazard_multiplier_m_rec: 3.2,
                systemic_crash_hazard_lambda_crash: 0.14 + credit_bump * 2.0,
                crash_state_probability_p_crash: 0.06,
                loss_given_default_lgd: 0.60,
                tranche_attachment_k1: 0.03,
                tranche_detachment_k2: 0.07,
                risk_free_rate_r: r,
            },
        );
        let gj_svi = ql_vol::evaluate_gatheral_jacquier_svi_raw_to_jump_wing_reparameterization(
            &ql_vol::GatheralJacquierSviJumpWingSpec {
                raw_svi_a: 0.012 + credit_bump * 0.5,
                raw_svi_b: 0.14,
                raw_svi_rho: -0.42,
                raw_svi_m: 0.05,
                raw_svi_sigma: 0.18 + credit_bump,
                evaluation_log_moneyness_k: -0.10,
                expiry_years_t: 0.75,
            },
        );
        let bk_xva = ql_xva::evaluate_burgard_kjaer_semi_replication_funding_strategy_xva(
            &ql_xva::BurgardKjaerSemiReplicationFundingStrategyXvaSpec {
                risk_free_derivative_value_v0: 100.0 + credit_bump * 100.0,
                negative_exposure_ene0: 45.0,
                counterparty_hazard_rate_lambda_c: 0.022 + credit_bump,
                bank_hazard_rate_lambda_b: 0.015,
                loss_given_default_lgd: 0.60,
                own_bond_recovery_r_b: 0.40,
                semi_replication_hedge_weight_alpha: 0.65,
                ois_discount_rate_r: r,
            },
        );
        let benth_hdd =
            ql_commodities::evaluate_benth_cartea_saltyte_benth_seasonal_temperature_hdd_cdd_swap(
                &ql_commodities::BenthSeasonalTemperatureHddCddSpec {
                    initial_temperature_t0: 12.5 - inflation * 10.0,
                    seasonal_mean_level_a0: 15.0,
                    seasonal_amplitude_b1: 11.0,
                    mean_reversion_speed_kappa: 6.5,
                    temperature_volatility_sigma: 8.0 + credit_bump * 20.0,
                    degree_day_threshold_k_base: 18.0,
                    market_price_of_weather_risk_theta: -0.18,
                    discount_rate_r: r,
                },
            );
        let bp_fx = ql_fx::evaluate_baba_packer_nagano_fx_swap_implied_dollar_funding_basis(
            &ql_fx::BabaPackerFxSwapDollarBasisSpec {
                spot_fx_s0: 1.085,
                outright_forward_fx_f: 1.096,
                usd_ois_rate_r_usd: r + 0.015,
                foreign_ois_rate_r_for: r,
                dealer_balance_sheet_leverage_cost_lambda_bs: 0.0022 + credit_bump * 0.5,
                foreign_dollar_funding_demand_pressure_d_usd: 0.0035,
                swap_tenor_years_t: 0.50,
            },
        );
        let cv_lb =
            ql_exotics::evaluate_conze_viswanathan_floating_and_fixed_lookback_extremum_option(
                &ql_exotics::ConzeViswanathanLookbackSpec {
                    spot_price_s0: spot,
                    observed_running_minimum_s_min: 88.0,
                    risk_free_rate_r: r,
                    dividend_yield_q: 0.015,
                    volatility_sigma: vol,
                    discrete_monitoring_steps_m: 52.0,
                    maturity_years_t: 1.0,
                },
            );
        let hl_bl = ql_portfolio::evaluate_litterman_he_institutional_equilibrium_posterior_black_litterman_tilt(
            &ql_portfolio::HeLittermanBlackLittermanTiltSpec {
                risk_aversion_delta: 2.5,
                market_cap_weight_asset1_w1: 0.60,
                asset1_volatility_sigma1: vol,
                asset2_volatility_sigma2: 0.22 + credit_bump,
                asset_correlation_rho: 0.35,
                prior_uncertainty_scalar_tau: 0.05,
                relative_view_outperformance_q: 0.035,
                view_error_volatility_sigma_omega: 0.04,
            },
        );
        let ghs_is = ql_portfolio::evaluate_glasserman_heidelberger_shahabuddin_importance_sampling_credit_var(
            &ql_portfolio::GlassermanImportanceSamplingVaRSpec {
                portfolio_linear_delta_exposure_a: 1.80 + credit_bump * 10.0,
                portfolio_quadratic_gamma_exposure_b: 0.22,
                idiosyncratic_noise_volatility_sigma_e: 0.65,
                tail_loss_threshold_x: 5.20,
                constant_loss_drift_c0: 0.15,
                confidence_level_alpha: 0.99,
                effective_sample_paths_n: 1000.0,
            },
        );
        let egj_vecm =
            ql_dml::evaluate_engle_granger_johansen_cointegration_error_correction_pairs_spread(
                &ql_dml::EngleGrangerJohansenCointegrationSpec {
                    current_log_price_y: 4.65 + credit_bump,
                    current_log_price_x: 4.52,
                    cointegration_hedge_ratio_beta: 1.04,
                    equilibrium_intercept_mu: -0.08,
                    error_correction_speed_alpha: 0.28,
                    innovation_volatility_sigma_eps: 0.032,
                    canonical_correlation_r_sq: 0.24,
                },
            );
        let gab_dis =
            ql_models::evaluate_gabaix_rare_disasters_time_varying_resilience_equity_premium(
                &ql_models::GabaixRareDisastersResilienceSpec {
                    subjective_discount_rate_rho: 0.032,
                    relative_risk_aversion_gamma: 4.0,
                    normal_consumption_growth_g_c: 0.020,
                    normal_consumption_volatility_sigma_c: 0.018,
                    disaster_probability_p: 0.017 + credit_bump,
                    consumption_recovery_fraction_b_c: 0.72,
                    dividend_recovery_fraction_b_d: 0.18,
                    resilience_mean_reversion_phi_h: 0.18,
                },
            );

        merc_fva.fva_adjusted_uncollateralized_swap_rate * 100.0
            + ing_cb.convertible_bond_fair_value * 0.1
            + hw_cdo.bespoke_tranche_par_spread * 100.0
            + gj_svi.svi_implied_volatility_sigma_k * 10.0
            + bk_xva.total_burgard_kjaer_xva_charge * 0.1
            + benth_hdd.discounted_hdd_swap_present_value * 0.01
            + bp_fx.balance_sheet_adjusted_dollar_funding_pv * 100.0
            + cv_lb.discrete_floating_lookback_call_price * 0.1
            + hl_bl.posterior_expected_return_asset1 * 100.0
            + ghs_is.saddlepoint_adjusted_quadratic_var
            + egj_vecm.cointegration_confidence_weighted_sharpe
            + gab_dis.gabaix_equilibrium_price_dividend_ratio
    }
    pub fn mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let ab_tb = ql_curves::evaluate_ametrano_bianchetti_dual_curve_tenor_basis_swap_bootstrap(
            &ql_curves::AmetranoBianchettiTenorBasisSpec {
                ois_discount_rate_r_ois: r,
                ois_slope_theta_ois: 0.004,
                ibor_3m_spread_over_ois_s3: 0.0018 + credit_bump * 0.5,
                tenor_basis_3m_vs_6m_spread_b36: 0.0014 + credit_bump * 0.5,
                tenor_basis_term_steepening_gamma: 0.0005,
                compounding_convexity_vol_sigma3: 0.0075,
                swap_maturity_years_t: 2.0,
            },
        );
        let fh_iopo =
            ql_fixed_income::evaluate_fabozzi_hayre_io_po_stripped_mbs_prepayment_convexity(
                &ql_fixed_income::FabozziHayreIoPoStripSpec {
                    pool_principal_balance: 100.0,
                    weighted_average_coupon_wac: 0.065,
                    net_passthrough_coupon_c_net: 0.060,
                    current_mortgage_market_rate_r_mkt: r + 0.015,
                    baseline_cpr_min: 0.06,
                    peak_refi_cpr_max: 0.42,
                    refi_s_curve_steepness_beta: 35.0,
                    option_adjusted_spread_oas: 0.008 + credit_bump,
                },
            );
        let hb_wwr = ql_credit::evaluate_hull_white_brigo_cds_counterparty_wrong_way_spread(
            &ql_credit::HullBrigoWwrCdsSpec {
                reference_entity_hazard_lambda_ref: 0.024 + credit_bump,
                protection_seller_hazard_lambda_cp: 0.016 + credit_bump * 0.5,
                default_correlation_rho_wwr: 0.45,
                reference_recovery_rate_r_ref: 0.40,
                counterparty_recovery_rate_r_cp: 0.35,
                cds_spread_volatility_sigma_s: 0.42,
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let brs_vv = ql_vol::evaluate_bossens_rayee_skovmand_vanna_volga_barrier_survival_overhedge(
            &ql_vol::BossensVannaVolgaBarrierSpec {
                spot_fx_s0: 1.10,
                strike_k: 1.12,
                down_and_out_barrier_h: 1.02,
                atm_volatility_sigma_atm: 0.10 + credit_bump,
                risk_reversal_25d_rr: -0.018,
                butterfly_25d_bf: 0.0045,
                domestic_foreign_rate_diff_mu: r - 0.02,
                maturity_years_t: 1.0,
            },
        );
        let greg_mva = ql_xva::evaluate_gregory_isda_simm_dynamic_initial_margin_mva(
            &ql_xva::GregoryIsdaSimmMvaSpec {
                trade_notional_dv01_w0: 120.0 + credit_bump * 100.0,
                simm_calibrated_volatility_sigma_im: 0.011,
                margin_period_of_risk_days_mpor: 10.0,
                im_funding_spread_s_im: 0.0085 + credit_bump,
                counterparty_hazard_rate_lambda_c: 0.018 + credit_bump,
                bank_hazard_rate_lambda_b: 0.012,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let cd_cds =
            ql_commodities::evaluate_carmona_durrleman_clean_dark_spread_carbon_power_option(
                &ql_commodities::CarmonaCleanDarkSpreadSpec {
                    electricity_forward_f_power: 95.0 + inflation * 50.0,
                    coal_forward_f_coal: 110.0,
                    carbon_allowance_forward_f_eua: 68.0,
                    coal_heat_rate_h_coal: 0.38,
                    carbon_emission_intensity_e_co2: 0.92,
                    fixed_om_strike_k: 4.50,
                    power_volatility_sigma_p: 0.38 + credit_bump,
                    fuel_carbon_basket_volatility_sigma_fc: 0.30,
                    power_fuel_correlation_rho: 0.72,
                },
            );
        let rh_qto = ql_fx::evaluate_reiner_haug_foreign_equity_domestic_strike_and_quanto_call(
            &ql_fx::ReinerHaugQuantoEquityCallSpec {
                foreign_equity_spot_s0: spot,
                fx_spot_e0: 1.10,
                foreign_strike_k_for: 102.0,
                domestic_rate_rd: r,
                foreign_rate_rf: 0.025,
                equity_dividend_yield_q: 0.015,
                equity_volatility_sigma_s: vol,
                fx_volatility_sigma_e: 0.12 + credit_bump,
                equity_fx_correlation_rho: -0.35,
            },
        );
        let hr_soft = ql_exotics::evaluate_hart_ross_haug_soft_proportional_band_barrier_option(
            &ql_exotics::HartRossSoftBarrierSpec {
                spot_price_s0: spot,
                strike_price_k: 100.0,
                lower_soft_barrier_l: 80.0,
                upper_soft_barrier_u: 90.0,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let ps_alp = ql_portfolio::evaluate_pastor_stambaugh_bayesian_mispricing_alpha_shrinkage(
            &ql_portfolio::PastorStambaughAlphaShrinkageSpec {
                ols_sample_alpha_hat: 0.042,
                prior_mispricing_uncertainty_sigma_alpha: 0.020,
                residual_idiosyncratic_volatility_sigma_eps: vol,
                sample_history_years_t_obs: 5.0,
                benchmark_factor_sharpe_ratio_sr_f: 0.45,
                asset_factor_beta: 1.10,
                benchmark_expected_excess_return_mu_f: 0.065 - credit_bump,
            },
        );
        let ow_lob =
            ql_portfolio::evaluate_obizhaeva_wang_transient_resilience_impulse_block_execution(
                &ql_portfolio::ObizhaevaWangImpulseExecutionSpec {
                    order_size_shares_x0: 100.0,
                    permanent_price_impact_lambda: 0.0015 + credit_bump * 0.02,
                    transient_lob_impact_kappa: 0.0045,
                    lob_resilience_speed_rho: 2.20,
                    execution_horizon_t: 1.0,
                    asset_volatility_sigma: vol,
                    inventory_risk_aversion_gamma: 0.02,
                },
            );
        let dmw_hln = ql_dml::evaluate_diebold_mariano_west_hln_forecast_comparison_test(
            &ql_dml::DieboldMarianoWestForecastSpec {
                mean_loss_differential_d_bar: 0.028 + credit_bump,
                lag0_loss_diff_variance_gamma0: 0.045,
                lag1_loss_diff_autocovariance_gamma1: 0.014,
                forecast_horizon_steps_h: 3.0,
                out_of_sample_observations_p: 120.0,
                in_sample_estimation_observations_r: 240.0,
                parameter_estimation_variance_v_theta: 0.008,
            },
        );
        let ps_nom = ql_models::evaluate_piazzesi_schneider_stagflation_nominal_bond_term_premium(
            &ql_models::PiazzesiSchneiderNominalBondSpec {
                real_short_rate_r_real: r * 0.5,
                expected_inflation_pi0: inflation,
                relative_risk_aversion_gamma: 8.5,
                real_growth_persistence_phi_x: 0.92,
                inflation_persistence_phi_pi: 0.94,
                real_growth_volatility_sigma_x: 0.006,
                inflation_volatility_sigma_pi: 0.008 + credit_bump * 0.2,
                growth_inflation_correlation_rho_x_pi: -0.55,
            },
        );

        ab_tb.bootstrapped_6m_ibor_par_swap_rate * 100.0
            + fh_iopo.io_strip_present_value * 0.5
            + hb_wwr.wrong_way_adjusted_cds_par_spread * 100.0
            + brs_vv.vanna_volga_barrier_call_price * 100.0
            + greg_mva.lifetime_simm_mva_charge * 10.0
            + cd_cds.clean_dark_spread_call_option_value * 0.5
            + rh_qto.fixed_fx_quanto_call_present_value * 0.5
            + hr_soft.soft_down_and_out_call_price * 0.5
            + ps_alp.posterior_predictive_information_ratio * 10.0
            + ow_lob.risk_adjusted_total_execution_cost * 0.1
            + dmw_hln.hln_adjusted_dmw_t_statistic
            + ps_nom.equilibrium_long_nominal_bond_yield * 100.0
    }
    pub fn ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let sw_ufr = ql_curves::evaluate_smith_wilson_eiopa_ultimate_forward_rate_extrapolation(
            &ql_curves::SmithWilsonUfrExtrapolationSpec {
                intermediate_pillar_years_u1: 10.0,
                last_liquid_point_years_u2: 20.0,
                extrapolated_maturity_years_t_ext: 50.0,
                intermediate_zero_rate_r1: r * 0.9,
                last_liquid_point_zero_rate_r2: r,
                ultimate_forward_rate_omega: 0.036 + inflation * 0.1,
                convergence_speed_alpha: 0.12,
            },
        );
        let kwf_cb =
            ql_fixed_income::evaluate_kalotay_williams_fabozzi_callable_bond_refunding_efficiency(
                &ql_fixed_income::KalotayWilliamsFabozziCallableSpec {
                    initial_short_rate_r0: r,
                    short_rate_drift_mu: 0.002,
                    lognormal_short_rate_volatility_sigma: 0.22,
                    option_adjusted_spread_oas: 0.010 + credit_bump,
                    coupon_rate_c: 0.058,
                    call_strike_price_k: 100.2,
                    refunding_transaction_cost: 0.35,
                    step_years_dt: 1.0,
                },
            );
        let sc_cdo = ql_credit::evaluate_schonbucher_archimedean_clayton_copula_cdo_tranche_loss(
            &ql_credit::SchonbucherClaytonCopulaCdoSpec {
                pool_hazard_rate_lambda: 0.022 + credit_bump,
                recovery_rate_r: 0.40,
                clayton_dependence_theta: 0.65,
                attachment_point_ka: 0.03,
                detachment_point_kb: 0.07,
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let bue_vc = ql_vol::evaluate_buehler_consistent_double_exponential_variance_curve(
            &ql_vol::BuehlerVarianceCurveSpec {
                short_variance_xi0: vol * vol,
                long_variance_xi_inf: 0.048 + credit_bump,
                hump_variance_amplitude_xi_hump: 0.015,
                fast_mean_reversion_kappa1: 2.4,
                slow_mean_reversion_kappa2: 0.55,
                fast_factor_vol_of_vol_omega1: 1.35,
                slow_factor_vol_of_vol_omega2: 0.45,
                factor_correlation_rho12: 0.35,
                maturity_years_t: 1.5,
            },
        );
        let pr_cva = ql_xva::evaluate_pykhtin_rosen_euler_marginal_cva_netting_allocation(
            &ql_xva::PykhtinRosenMarginalCvaSpec {
                subportfolio_1_expected_mtm_mu1: 12.0,
                subportfolio_2_expected_mtm_mu2: -4.5,
                subportfolio_1_volatility_sigma1: 28.0,
                subportfolio_2_volatility_sigma2: 22.0,
                mtm_correlation_rho12: -0.35,
                counterparty_hazard_lambda_c: 0.024 + credit_bump,
                counterparty_recovery_rate_rc: 0.40,
                risk_free_rate_r: r,
            },
        );
        let hr_fut =
            ql_commodities::evaluate_hilliard_reis_stochastic_convenience_yield_jump_futures_option(
                &ql_commodities::HilliardReisJumpFuturesSpec {
                    futures_price_f0: 78.0 + inflation * 20.0,
                    strike_price_k: 80.0,
                    spot_volatility_sigma_s: 0.32 + credit_bump,
                    convenience_yield_volatility_sigma_delta: 0.24,
                    convenience_yield_mean_reversion_kappa: 1.60,
                    spot_convenience_correlation_rho: 0.65,
                    jump_intensity_lambda_j: 0.85,
                    jump_volatility_v_j: 0.18,
                    risk_free_rate_r: r,
                },
            );
        let gm_fx =
            ql_fx::evaluate_gabaix_maggiori_imperfect_financial_markets_fx_carry_risk_premium(
                &ql_fx::GabaixMaggioriFxCarrySpec {
                    home_gross_interest_rate_r: 1.0 + r,
                    foreign_gross_interest_rate_r_star: 1.015,
                    financier_friction_gamma: 0.85 + credit_bump * 5.0,
                    fx_volatility_sigma: 0.14,
                    capital_flow_demand_f0: 0.22,
                    trade_flow_elasticity_iota0: 1.00,
                    period_1_trade_elasticity_iota1: 1.05,
                },
            );
        let cv_shout = ql_exotics::evaluate_cheuk_vorst_haug_strike_reset_shout_call_option(
            &ql_exotics::CheukVorstShoutResetSpec {
                spot_s0: spot,
                initial_strike_k0: 102.0,
                reset_monitoring_time_t1: 0.45,
                maturity_years_t: 1.0,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
            },
        );
        let gz_kelly =
            ql_portfolio::evaluate_maclean_thorp_ziemba_grossman_zhou_fractional_kelly_drawdown(
                &ql_portfolio::GrossmanZhouDrawdownKellySpec {
                    expected_return_mu: 0.115,
                    risk_free_rate_r: r,
                    asset_volatility_sigma: vol,
                    max_drawdown_floor_delta: 0.25,
                    cushion_risk_aversion_gamma: 2.0,
                    current_drawdown_from_peak: 0.06 + credit_bump,
                },
            );
        let gls_raroc = ql_portfolio::evaluate_gourieroux_laurent_scaillet_tasche_cornish_fisher_euler_var_raroc(
            &ql_portfolio::GlsTascheCornishFisherRarocSpec {
                portfolio_expected_excess_return_mu_p: 0.082,
                portfolio_volatility_sigma_p: 0.145 + credit_bump,
                portfolio_skewness_s_p: -0.65,
                portfolio_excess_kurtosis_k_p: 2.40,
                normal_tail_quantile_z_alpha: 2.326_347_874,
                subportfolio_expected_excess_return_mu_i: 0.105,
                subportfolio_volatility_sigma_i: vol,
                subportfolio_portfolio_correlation_rho: 0.72,
            },
        );
        let ap_wald =
            ql_dml::evaluate_andrews_ploberger_hansen_sup_exp_ave_wald_structural_break_test(
                &ql_dml::AndrewsPlobergerStructuralBreakSpec {
                    sample_size_t_obs: 120.0,
                    pre_break_slope_beta1: 0.45,
                    post_break_slope_beta2: 0.72 + credit_bump * 2.0,
                    break_profile_dispersion_delta: 0.14,
                    residual_variance_sigma2: 0.55,
                    regressor_second_moment_q_xx: 1.15,
                    trimming_fraction_pi0: 0.15,
                },
            );
        let jq_fin = ql_models::evaluate_jermann_quadrini_financial_shocks_enforcement_constraint_credit_cycle(
            &ql_models::JermannQuadriniFinancialShockSpec {
                capital_stock_k: 10.0,
                working_capital_output_y: 2.4,
                financial_enforcement_xi: 0.42 - credit_bump,
                risk_free_rate_r: r,
                corporate_tax_shield_tau: 0.28,
                household_discount_factor_beta: 0.960,
                equity_payout_cost_phi: 0.18,
                equity_payout_gap_d: -0.15,
            },
        );

        sw_ufr.extrapolated_zero_yield * 100.0
            + kwf_cb.refunding_efficiency_ratio * 10.0
            + sc_cdo.fair_tranche_par_spread * 100.0
            + bue_vc.fair_variance_swap_volatility_strike * 10.0
            + pr_cva.total_netting_set_cva
            + hr_fut.jump_diffusion_futures_call_price * 0.5
            + gm_fx.expected_fx_carry_excess_return * 100.0
            + cv_shout.strike_reset_call_value * 0.5
            + gz_kelly.expected_geometric_growth_rate * 100.0
            + gls_raroc.subportfolio_euler_raroc * 10.0
            + ap_wald.exponential_wald_statistic_exp_w
            + jq_fin.endogenous_financial_labor_wedge * 100.0
    }
    pub fn smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let mc_tax = ql_curves::evaluate_mcculloch_litzenberger_tax_clientele_discount_curve(
            &ql_curves::McCullochTaxClienteleCurveSpec {
                low_coupon_rate_c_low: 0.025,
                high_coupon_rate_c_high: 0.060,
                low_coupon_bond_price_p_low: 0.985,
                high_coupon_bond_price_p_high: 1.010,
                marginal_income_tax_rate_tau_inc: 0.28 + inflation * 0.2,
                capital_gains_tax_rate_tau_cg: 0.15,
                municipal_credit_liquidity_spread: 0.0035 + credit_bump,
            },
        );
        let ab_cb = ql_fixed_income::evaluate_andersen_buffum_equity_linked_hazard_convertible_bond(
            &ql_fixed_income::AndersenBuffumConvertibleBondSpec {
                spot_price_s0: spot,
                conversion_ratio: 1.05,
                baseline_hazard_lambda0: 0.028 + credit_bump,
                hazard_power_elasticity_p: 1.25,
                bond_recovery_rate_r: 0.40,
                coupon_rate_c: 0.045,
                risk_free_rate_r: r,
                stock_volatility_sigma: vol,
            },
        );
        let ok_rl = ql_credit::evaluate_okane_brigo_recovery_lock_and_digital_default_swap(
            &ql_credit::OkaneBrigoRecoveryLockCdsSpec {
                market_cds_spread_s_mkt: 0.0180 + credit_bump,
                conventional_quoting_recovery_r_conv: 0.40,
                expected_physical_recovery_r_bar: 0.34,
                recovery_lock_strike_r_lock: 0.38,
                hazard_recovery_correlation_rho: -0.45,
                recovery_volatility_sigma_r: 0.18,
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let dk_tree = ql_vol::evaluate_derman_kani_implied_tree_arrow_debreu_local_volatility(
            &ql_vol::DermanKaniImpliedTreeSpec {
                central_node_spot_si: spot,
                arrow_debreu_state_price_lambda: 0.96,
                atm_implied_volatility_sigma: vol,
                smile_skew_slope_beta: -0.15 - credit_bump,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                step_years_dt: 0.25,
            },
        );
        let kg_lva = ql_xva::evaluate_kenyon_green_asymmetric_csa_collateral_liquidity_lva(
            &ql_xva::KenyonGreenCsaLvaSpec {
                expected_portfolio_mtm_mu: 6.0,
                portfolio_mtm_volatility_sigma: 24.0 + credit_bump * 50.0,
                csa_receive_threshold_h_rec: 5.0,
                csa_post_threshold_h_post: 2.0,
                posted_collateral_liquidity_spread: 0.0045 + credit_bump,
                received_collateral_remuneration_drag: 0.0020,
                ois_discount_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let bdj_gas = ql_commodities::evaluate_boogert_de_jong_gas_storage_ratchet_injection_withdrawal_spread(
            &ql_commodities::BoogertDeJongGasRatchetSpec {
                summer_injection_forward_f_sum: 2.85,
                winter_withdrawal_forward_f_win: 4.10 + inflation * 5.0,
                working_gas_fill_fraction: 0.45,
                ratchet_pressure_attenuation_alpha: 0.35,
                variable_injection_withdrawal_cost: 0.12,
                calendar_spread_normal_volatility_sigma: 0.85 + credit_bump * 2.0,
                risk_free_rate_r: r,
            },
        );
        let gk_baw = ql_fx::evaluate_garman_kohlhagen_biger_hull_fx_early_exercise_premium(
            &ql_fx::GarmanKohlhagenBawFxSpec {
                spot_fx_s0: 1.10,
                strike_k: 1.08,
                domestic_rate_rd: r,
                foreign_rate_rf: 0.055 + credit_bump,
                fx_volatility_sigma: 0.12,
                maturity_years_t: 1.0,
            },
        );
        let hh_bin = ql_exotics::evaluate_hui_haug_one_touch_and_no_touch_binary_barrier(
            &ql_exotics::HuiHaugBinaryBarrierSpec {
                spot_price_s0: spot,
                down_barrier_h: spot * 0.88,
                cash_rebate_k: 10.0,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let kr_hj = ql_portfolio::evaluate_kan_robotti_misspecification_robust_hj_factor_distance(
            &ql_portfolio::KanRobottiHjDistanceSpec {
                traded_factor_mean_mu1: 0.068,
                nontraded_factor_price_lambda2: 0.045,
                test_asset_1_pricing_error_alpha1: 0.018 + credit_bump,
                test_asset_2_pricing_error_alpha2: -0.014,
                test_asset_1_variance_v1: 0.036,
                test_asset_2_variance_v2: 0.044,
                test_asset_correlation_rho12: 0.30,
                sample_observations_t_obs: 180.0,
            },
        );
        let cj_of = ql_portfolio::evaluate_cartea_jaimungal_order_flow_alpha_optimal_liquidation(
            &ql_portfolio::CarteaJaimungalOrderFlowLiquidationSpec {
                initial_inventory_q0: 100.0,
                order_flow_signal_mu0: 0.35 - credit_bump * 2.0,
                signal_mean_reversion_kappa: 2.50,
                temporary_impact_k: 0.0040,
                running_inventory_penalty_phi: 0.012,
                terminal_liquidation_penalty_a: 0.080,
                horizon_t: 1.0,
            },
        );
        let asj_jmp = ql_dml::evaluate_ait_sahalia_jacod_power_variation_jump_activity_test(
            &ql_dml::AitSahaliaJacodJumpTestSpec {
                fine_grid_power_variation_b_fine: 0.0085 + credit_bump * 0.02,
                coarse_grid_power_variation_b_coarse: 0.0092,
                power_exponent_p: 4.0,
                subsampling_multiplier_k: 2.0,
                multipower_quarticity_scale_a2p: 0.00018,
                high_frequency_sample_size_n: 390.0,
            },
        );
        let wac_dis = ql_models::evaluate_wachter_time_varying_disaster_intensity_equity_volatility(
            &ql_models::WachterDisasterVolatilitySpec {
                current_disaster_intensity_lambda: 0.0355 + credit_bump,
                long_run_disaster_intensity_lambda_bar: 0.0355,
                intensity_mean_reversion_kappa: 0.080,
                intensity_volatility_sigma_lambda: 0.067,
                diffusive_consumption_volatility_sigma_c: 0.020,
                relative_risk_aversion_gamma: 3.0,
                disaster_contraction_fraction_b: 0.22,
                dividend_leverage_phi: 2.6,
            },
        );

        mc_tax.municipal_to_treasury_yield_ratio * 10.0
            + ab_cb.convertible_bond_value * 0.1
            + ok_rl.recovery_adjusted_digital_protection_pv * 100.0
            + dk_tree.implied_node_local_volatility * 10.0
            + kg_lva.total_asymmetric_csa_lva_charge * 10.0
            + bdj_gas.ratchet_adjusted_storage_option_value * 10.0
            + gk_baw.american_fx_call_value * 100.0
            + hh_bin.cash_at_hit_one_touch_value
            + kr_hj.kan_robotti_hj_t_statistic * 10.0
            + cj_of.optimal_initial_trading_speed * 0.1
            + asj_jmp.standardized_asj_jump_z_statistic
            + wac_dis.equilibrium_equity_risk_premium * 100.0
    }
    pub fn mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let hw_mc = ql_curves::evaluate_hagan_west_monotone_convex_forward_interpolation(
            &ql_curves::HaganWestMonotoneConvexSpec {
                interval_start_years_t0: 2.0,
                interval_end_years_t1: 5.0,
                target_maturity_years_t: 3.4,
                start_zero_rate_r0: r,
                previous_discrete_forward_f_prev: r + 0.001,
                current_discrete_forward_f_curr: r + 0.007,
                next_discrete_forward_f_next: r + 0.013,
                node_blending_weight_w: 0.50,
            },
        );
        let fab_oas =
            ql_fixed_income::evaluate_fabozzi_tuckman_effective_duration_and_convexity_oas(
                &ql_fixed_income::FabozziEffectiveDurationOasSpec {
                    benchmark_yield_y0: r,
                    option_adjusted_spread_oas: 0.0085 + credit_bump,
                    coupon_rate_c: 0.056,
                    call_price_k: 100.5,
                    bond_price_volatility_sigma: 0.065,
                    parallel_yield_shift_dy: 0.0025,
                    maturity_years_t: 5.0,
                },
            );
        let fa_lgd = ql_credit::evaluate_frye_altman_procyclical_downturn_lgd_credit_loss(
            &ql_credit::FryeAltmanDownturnLgdSpec {
                default_threshold_z_pd: -2.0537 + credit_bump * 5.0,
                asset_factor_correlation_rho: 0.18,
                expected_recovery_rate_mu_r: 0.45,
                recovery_volatility_sigma_r: 0.20,
                recovery_factor_loading_q_r: 0.55,
                adverse_macro_quantile_z_alpha: 2.3263,
            },
        );
        let gj_svi = ql_vol::evaluate_gatheral_jacquier_svi_calendar_and_butterfly_arbitrage(
            &ql_vol::GatheralJacquierSviSurfaceArbitrageSpec {
                log_moneyness_k: -0.10,
                svi_level_a: vol * vol * 0.6,
                svi_slope_b: 0.14 + credit_bump,
                svi_correlation_rho: -0.45,
                svi_center_m: 0.02,
                svi_smoothness_sigma: 0.18,
                maturity_years_t1: 1.0,
                calendar_variance_growth_dw_dt: 0.048,
            },
        );
        let gb_bcva = ql_xva::evaluate_gregory_brigo_first_to_default_bilateral_cva_dva(
            &ql_xva::GregoryBrigoBilateralCvaDvaSpec {
                expected_positive_exposure_epe: 18.0,
                expected_negative_exposure_ene: 11.0,
                counterparty_hazard_lambda_c: 0.026 + credit_bump,
                bank_own_hazard_lambda_b: 0.015,
                common_shock_hazard_lambda_cb: 0.004,
                counterparty_lgd_c: 0.60,
                bank_own_lgd_b: 0.60,
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let bs_mine =
            ql_commodities::evaluate_brennan_schwartz_natural_resource_mine_hysteresis_option(
                &ql_commodities::BrennanSchwartzMineHysteresisSpec {
                    commodity_spot_price_p0: 4.10 + inflation * 5.0,
                    variable_extraction_cost_c: 2.40,
                    reopening_investment_cost_i: 8.50,
                    suspension_closure_cost_k: 1.80,
                    convenience_yield_delta: 0.045,
                    risk_free_rate_r: r,
                    commodity_volatility_sigma: 0.26 + credit_bump,
                },
            );
        let wys_fad = ql_fx::evaluate_wystup_fx_corridor_fader_option_and_forward(
            &ql_fx::WystupFxFaderOptionSpec {
                spot_fx_s0: 1.09,
                strike_k: 1.07,
                corridor_lower_bound_l: 1.03,
                corridor_upper_bound_u: 1.16,
                domestic_rate_rd: r,
                foreign_rate_rf: 0.024,
                fx_volatility_sigma: 0.11 + credit_bump,
                maturity_years_t: 1.0,
            },
        );
        let bd_cap = ql_exotics::evaluate_broadie_detemple_american_capped_call_option(
            &ql_exotics::BroadieDetempleAmericanCappedCallSpec {
                spot_price_s0: spot,
                strike_k: spot * 0.96,
                cap_level_l: spot * 1.18,
                risk_free_rate_r: r,
                dividend_yield_q: 0.025,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let meu_ep =
            ql_portfolio::evaluate_meucci_satchell_scowcroft_entropy_pooling_view_allocation(
                &ql_portfolio::MeucciEntropyPoolingAllocationSpec {
                    prior_mean_asset1_mu1: 0.065,
                    prior_mean_asset2_mu2: 0.050,
                    prior_variance_asset1_v1: 0.040,
                    prior_variance_asset2_v2: 0.030,
                    prior_correlation_rho12: 0.35,
                    spread_view_target_v: 0.032 + credit_bump,
                    spread_view_variance_tau2: 0.022,
                    view_confidence_c: 0.65,
                },
            );
        let cds_rob = ql_portfolio::evaluate_cont_deguest_scandolo_tail_risk_influence_sensitivity(
            &ql_portfolio::ContDeguestScandoloRobustRiskSpec {
                loss_location_mu: 0.012,
                tail_scale_s: 0.025 + credit_bump,
                tail_index_xi: 0.24,
                tail_probability_p: 0.025,
                contamination_epsilon: 0.002,
                stress_outlier_loss_z: 0.28,
            },
        );
        let pss_ardl = ql_dml::evaluate_pesaran_shin_smith_ardl_bounds_cointegration_test(
            &ql_dml::PesaranShinSmithArdlBoundsSpec {
                lagged_dependent_coef_pi_yy: -0.24,
                lagged_regressor_coef_pi_yx: 0.36 + credit_bump,
                lagged_y_variance_s_yy: 1.45,
                lagged_x_variance_s_xx: 1.80,
                lagged_yx_covariance_s_yx: 0.95,
                residual_variance_sigma2: 0.085,
                sample_size_t_obs: 120.0,
            },
        );
        let cee_nk = ql_models::evaluate_christiano_eichenbaum_evans_calvo_investment_nk_impulse(
            &ql_models::ChristianoEichenbaumEvansNkSpec {
                calvo_price_stickiness_xi_p: 0.68,
                inflation_indexation_iota_p: 0.45,
                discount_factor_beta: 0.99,
                habit_persistence_h: 0.70,
                investment_adjustment_curvature_kappa_i: 4.50,
                steady_state_investment_share: 0.22,
                real_interest_rate_shock: 0.010 + credit_bump,
            },
        );

        hw_mc.interpolated_zero_rate * 100.0
            + fab_oas.effective_duration
            + fa_lgd.procyclical_downturn_stressed_loss * 100.0
            + gj_svi.dupire_local_volatility * 10.0
            + gb_bcva.net_bilateral_cva_charge
            + bs_mine.operating_mine_total_value * 0.1
            + wys_fad.fade_in_corridor_call_value * 100.0
            + bd_cap.american_capped_call_value * 0.5
            + meu_ep.posterior_squared_sharpe_ratio * 10.0
            + cds_rob.contaminated_median_to_expected_shortfall_ratio * 10.0
            + pss_ardl.long_run_multiplier_standard_error * 10.0
            + cee_nk.contractionary_inflation_impulse_response * 100.0
    }
    pub fn hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let mr_basis = ql_curves::evaluate_mercurio_russo_ibor_ois_credit_liquidity_basis(
            &ql_curves::MercurioRussoIborOisCreditLiquiditySpec {
                ois_rate_r: r,
                interbank_credit_hazard_lambda: 0.0045 + credit_bump * 0.5,
                interbank_loss_given_default: 0.60,
                term_liquidity_spread_l: 0.0018 + credit_bump * 0.2,
                credit_liquidity_covariance_gamma: 0.0008,
                accrual_tenor_tau: 0.50,
            },
        );
        let jy_yoy =
            ql_fixed_income::evaluate_jarrow_yildirim_yoy_inflation_swap_convexity_adjustment(
                &ql_fixed_income::JarrowYildirimYoyInflationConvexitySpec {
                    nominal_forward_rate_f_n: r + 0.006,
                    real_forward_rate_f_r: r - inflation,
                    nominal_rate_volatility_sigma_n: 0.0085,
                    real_rate_volatility_sigma_r: 0.0065,
                    cpi_index_volatility_sigma_i: 0.012 + credit_bump,
                    nominal_real_correlation_rho_nr: 0.55,
                    real_cpi_correlation_rho_ri: -0.25,
                    period_start_years_t_prev: 4.0,
                },
            );
        let dl_cds = ql_credit::evaluate_duffie_lando_das_sundaram_imperfect_information_cds_spread(
            &ql_credit::DuffieLandoDasSundaramHybridCdsSpec {
                reported_log_distance_to_default: 0.55 - credit_bump * 4.0,
                asset_volatility_sigma_v: 0.25,
                accounting_noise_std_a: 0.18,
                surprise_jump_hazard_lambda0: 0.032 + credit_bump,
                jump_leverage_sensitivity_eta: 1.40,
                recovery_rate_r: 0.40,
                maturity_years_t: 1.0,
            },
        );
        let bfg_skew = ql_vol::evaluate_bayer_friz_gatheral_rough_bergomi_power_law_skew(
            &ql_vol::BayerFrizGatheralRoughSkewSpec {
                initial_forward_variance_xi0: vol * vol,
                hurst_exponent_h: 0.12,
                vol_of_vol_eta: 1.85 + credit_bump * 5.0,
                spot_vol_correlation_rho: -0.75,
                log_moneyness_k: -0.08,
                maturity_years_t: 0.25,
            },
        );
        let aps_mpor = ql_xva::evaluate_andersen_pykhtin_sokol_mpor_dispute_spike_cva(
            &ql_xva::AndersenPykhtinSokolMporDisputeSpec {
                portfolio_mtm_volatility_sigma: 42.0,
                base_mpor_years_delta: 0.040,
                unpaid_cashflow_spike_c: 3.20,
                variation_margin_receipt_fraction_alpha: 0.25,
                initial_margin_cushion_i0: 4.50,
                counterparty_hazard_lambda: 0.024 + credit_bump,
                counterparty_lgd: 0.60,
                ois_discount_rate_r: r,
            },
        );
        let bs_crack =
            ql_commodities::evaluate_bjerksund_stensland_three_asset_321_refinery_crack_option(
                &ql_commodities::BjerksundStenslandRefineryCrack321Spec {
                    crude_futures_price_f_c: 76.0 + inflation * 50.0,
                    gasoline_futures_price_f_g: 94.0 + inflation * 60.0,
                    distillate_futures_price_f_h: 98.0 + inflation * 65.0,
                    refining_strike_k: 16.0,
                    refined_basket_volatility_sigma_p: 0.29 + credit_bump,
                    crude_volatility_sigma_c: 0.27,
                    basket_crude_correlation_rho: 0.86,
                    risk_free_rate_r: r,
                },
            );
        let cc_pcs = ql_fx::evaluate_carr_chou_put_call_symmetry_static_barrier_replication(
            &ql_fx::CarrChouSymmetryBarrierReplicationSpec {
                spot_fx_s0: 1.12,
                strike_k: 1.10,
                down_barrier_h: 1.02,
                domestic_rate_rd: r,
                foreign_rate_rf: 0.028,
                call_strike_volatility_sigma_k: 0.115 + credit_bump,
                reflected_put_volatility_sigma_star: 0.132 + credit_bump,
                maturity_years_t: 1.0,
            },
        );
        let cjy_par = ql_exotics::evaluate_chesney_jeanblanc_yor_parisian_excursion_barrier_call(
            &ql_exotics::ChesneyJeanblancYorParisianBarrierSpec {
                spot_price_s0: spot,
                strike_k: spot * 0.95,
                down_barrier_h: spot * 0.87,
                excursion_window_years_d: 0.06,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
            },
        );
        let dmb_norm =
            ql_portfolio::evaluate_demiguel_brodie_norm_constrained_minimum_variance_shrinkage(
                &ql_portfolio::DeMiguelBrodieNormConstrainedPortfolioSpec {
                    asset_1_variance_v1: vol * vol,
                    asset_2_variance_v2: 0.028,
                    asset_correlation_rho12: 0.25,
                    ridge_norm_penalty_lambda2: 0.018 + credit_bump,
                    gross_exposure_penalty_lambda1: 0.008,
                    asset_1_expected_return_mu1: 0.075,
                    asset_2_expected_return_mu2: 0.052,
                },
            );
        let glft_mm = ql_portfolio::evaluate_gueant_lehalle_tapia_inventory_market_making_quotes(
            &ql_portfolio::GueantLehalleTapiaMarketMakingSpec {
                current_inventory_q: 3.0,
                mid_price_volatility_sigma: 0.45 + credit_bump * 5.0,
                risk_aversion_gamma: 0.12,
                order_intensity_decay_kappa: 1.60,
                baseline_arrival_intensity_a: 18.0,
                remaining_horizon_tau: 0.50,
            },
        );
        let hw_spa = ql_dml::evaluate_hansen_white_superior_predictive_ability_test(
            &ql_dml::HansenWhiteSuperiorPredictiveAbilitySpec {
                model_1_mean_loss_reduction_d1: 0.045 + credit_bump,
                model_2_mean_loss_reduction_d2: -0.024,
                model_1_hac_variance_w1: 0.065,
                model_2_hac_variance_w2: 0.080,
                cross_model_correlation_rho12: 0.35,
                out_of_sample_observations_n: 150.0,
            },
        );
        let bs_par = ql_models::evaluate_brunnermeier_sannikov_intermediary_volatility_paradox(
            &ql_models::BrunnermeierSannikovVolatilityParadoxSpec {
                intermediary_wealth_share_eta: 0.24 - credit_bump * 2.0,
                fundamental_capital_volatility_sigma: 0.025 + credit_bump * 0.5,
                expert_productivity_wedge_delta_a: 0.055,
                fire_sale_illiquidity_sensitivity_phi: 2.20,
                household_discount_rate_rho: 0.050,
                intermediary_leverage_ratio: 3.20,
            },
        );

        mr_basis.total_ibor_ois_additive_basis * 1000.0
            + jy_yoy.discounted_yoy_swaplet_value
            + dl_cds.fair_short_maturity_cds_spread * 100.0
            + bfg_skew.rough_bergomi_implied_volatility * 10.0
            + aps_mpor.aps_mpor_adjusted_cva
            + bs_crack.refinery_321_crack_call_value
            + cc_pcs.static_replicated_down_and_out_call * 100.0
            + cjy_par.parisian_down_and_out_call_value * 0.5
            + dmb_norm.regularized_portfolio_sharpe_ratio * 10.0
            + glft_mm.expected_instantaneous_spread_capture
            + hw_spa.hansen_spa_adjusted_p_value * 100.0
            + bs_par.volatility_paradox_amplification_ratio
    }
    pub fn mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let flt_lr = ql_curves::evaluate_filipovic_larsson_trolle_linear_rational_zlb_curve(
            &ql_curves::FilipovicLarssonTrolleLinearRationalCurveSpec {
                asymptotic_rate_alpha: r + 0.010,
                factor_1_state_x1: 0.35 + credit_bump * 2.0,
                factor_2_state_x2: 0.25,
                factor_1_mean_theta1: 0.30,
                factor_2_mean_theta2: 0.20,
                factor_1_speed_kappa1: 0.45,
                factor_2_speed_kappa2: 0.18,
                factor_1_volatility_sigma1: 0.14 + credit_bump,
                maturity_years_t: 5.0,
            },
        );
        let dk_repo =
            ql_fixed_income::evaluate_duffie_krishnamurthy_on_the_run_special_repo_premium(
                &ql_fixed_income::DuffieKrishnamurthySpecialRepoTreasurySpec {
                    general_collateral_repo_rate: r,
                    initial_repo_specialness_s0: 0.0110 + credit_bump,
                    specialness_decay_rate_kappa: 3.20,
                    short_squeeze_sensitivity_gamma: 0.0085,
                    short_to_float_ratio_q: 0.35,
                    remaining_auction_cycle_years: 0.25,
                    modified_duration_d_mod: 7.80,
                    dealer_intermediation_haircut_h: 0.18,
                },
            );
        let dsw_def = ql_credit::evaluate_duffie_saita_wang_stochastic_covariate_survival(
            &ql_credit::DuffieSaitaWangMultiPeriodDefaultSpec {
                baseline_log_hazard_beta0: -2.40 + credit_bump * 5.0,
                initial_distance_to_default_dd0: 2.60 - credit_bump * 5.0,
                long_run_distance_to_default_bar: 2.20,
                dd_mean_reversion_kappa: 0.38,
                dd_volatility_sigma: 0.45,
                dd_hazard_loading_beta_d: 0.85,
                horizon_years_t: 3.0,
            },
        );
        let cmn_log = ql_vol::evaluate_carr_madan_neuberger_log_contract_variance_replication(
            &ql_vol::CarrMadanNeubergerLogContractReplicationSpec {
                atm_implied_volatility_sigma: vol,
                smile_skew_slope_psi: -0.22 - credit_bump,
                smile_wing_curvature_kappa: 0.45,
                strike_grid_spacing_dk: 0.025,
                discrete_rebalancing_count_n: 252.0,
                maturity_years_t: 1.0,
            },
        );
        let hwbk_fva = ql_xva::evaluate_hull_white_burgard_kjaer_fva_dva2_separation(
            &ql_xva::HullWhiteBurgardKjaerFvaWindfallSpec {
                expected_positive_exposure_epe: 22.0,
                expected_negative_exposure_ene: 9.5,
                bank_default_hazard_lambda_b: 0.016 + credit_bump,
                bank_loss_given_default_lgd_b: 0.60,
                pure_liquidity_funding_spread_s_liq: 0.0045 + credit_bump * 0.5,
                counterparty_hazard_lambda_c: 0.020 + credit_bump,
                ois_risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let ccd_cy =
            ql_commodities::evaluate_casassus_collin_dufresne_level_dependent_convenience_yield(
                &ql_commodities::CasassusCollinDufresneConvenienceYieldSpec {
                    spot_price_s0: 82.0 + inflation * 50.0,
                    equilibrium_spot_anchor_s_bar: 70.0,
                    baseline_convenience_yield_delta0: 0.035,
                    spot_level_sensitivity_alpha_x: 0.42,
                    rate_sensitivity_alpha_r: 0.30,
                    risk_free_rate_r0: r,
                    spot_volatility_sigma_s: 0.28 + credit_bump,
                    market_price_of_risk_lambda_s: 0.35,
                    maturity_years_t: 1.5,
                },
            );
        let lrv_fx = ql_fx::evaluate_lustig_roussanov_verdelhan_fx_carry_factor_pricing(
            &ql_fx::LustigRoussanovVerdelhanFxCarryFactorSpec {
                interest_rate_differential_di: 0.038 + credit_bump,
                dollar_factor_beta_dol: 0.95,
                carry_slope_beta_hml: 0.65,
                global_vol_beta_vol: -0.42,
                dollar_factor_risk_price_lambda_dol: 0.015,
                carry_slope_risk_price_lambda_hml: 0.048,
                global_vol_risk_price_lambda_vol: -0.025,
                carry_factor_volatility_sigma_hml: 0.095 + credit_bump,
            },
        );
        let cl_corr = ql_exotics::evaluate_carr_lee_corridor_conditional_variance_swap(
            &ql_exotics::CarrLeeCorridorVarianceSwapSpec {
                forward_price_f0: spot,
                lower_corridor_barrier_b: spot * 0.94,
                atmf_volatility_sigma0: vol,
                skew_slope_psi: -0.18 - credit_bump,
                risk_free_rate_r: r,
                maturity_years_t: 1.0,
            },
        );
        let guw_mp = ql_portfolio::evaluate_garlappi_uppal_wang_multi_prior_ambiguity_portfolio(
            &ql_portfolio::GarlappiUppalWangMultiPriorPortfolioSpec {
                sample_mean_asset1_mu1: 0.085,
                sample_mean_asset2_mu2: 0.050,
                variance_asset1_v1: vol * vol,
                variance_asset2_v2: 0.026,
                correlation_rho12: 0.25,
                risk_aversion_gamma: 3.0,
                ambiguity_radius_epsilon: 0.42 + credit_bump * 5.0,
                estimation_sample_size_t: 60.0,
            },
        );
        let mf_evt = ql_portfolio::evaluate_mcneil_frey_garch_filtered_evt_pot_tail_risk(
            &ql_portfolio::McNeilFreyConditionalEvtGarchRiskSpec {
                conditional_loss_mean_mu: 0.0005,
                garch_omega: 0.000012,
                garch_arch_shock_term: 0.00018 + credit_bump * 0.01,
                garch_lagged_variance_term: 0.00022,
                standardized_pot_threshold_u: 1.65,
                threshold_exceedance_rate_pu: 0.060,
                gpd_shape_index_xi: 0.22 + credit_bump,
                gpd_scale_beta: 0.52,
            },
        );
        let pat_sjc = ql_dml::evaluate_patton_symmetrized_joe_clayton_tail_dependence(
            &ql_dml::PattonSymmetrizedJoeClaytonCopulaSpec {
                lower_tail_intercept_omega_l: -0.35 + credit_bump * 5.0,
                lower_tail_persistence_beta_l: 1.40,
                lagged_lower_tail_tau_l: 0.38,
                upper_tail_intercept_omega_u: -0.85,
                upper_tail_persistence_beta_u: 1.10,
                lagged_upper_tail_tau_u: 0.24,
                pit_discrepancy_forcing_delta: 0.18,
                forcing_sensitivity_alpha: -1.60,
            },
        );
        let pv_pol = ql_models::evaluate_pastor_veronesi_political_uncertainty_risk_premium(
            &ql_models::PastorVeronesiPoliticalUncertaintySpec {
                perceived_policy_impact_g_hat: -0.012 - credit_bump,
                prior_policy_uncertainty_sigma_g0: 0.035 + credit_bump,
                fundamental_macro_volatility_sigma_m: 0.14,
                political_cost_dispersion_sigma_pol: 0.085,
                investor_risk_aversion_gamma: 3.50,
                bayesian_learning_duration_tau: 2.0,
            },
        );

        flt_lr.state_dependent_yield_volatility * 10.0
            + dk_repo.annualized_specialness_pickup_rate
            + dsw_def.multi_period_default_probability * 10.0
            + cmn_log.discrete_replication_rmse * 1000.0
            + hwbk_fva.funding_benefit_adjustment_fba
            + ccd_cy.equilibrium_futures_price * 0.1
            + lrv_fx.carry_factor_sharpe_ratio * 10.0
            + cl_corr.corridor_variance_swap_pv
            + guw_mp.max_min_certainty_equivalent * 100.0
            + mf_evt.conditional_evt_garch_es99 * 100.0
            + pat_sjc.symmetrized_tail_asymmetry_wedge * 10.0
            + pv_pol.political_uncertainty_equity_risk_premium * 100.0
    }
    pub fn filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let lm_rfr = ql_curves::evaluate_lyashenko_mercurio_rfr_in_arrears_caplet_convexity(
            &ql_curves::LyashenkoMercurioRfrInArrearsCapletSpec {
                forward_rfr_rate_f0: r + 0.010,
                strike_rate_k: 0.040,
                fixing_start_years_t_prev: 0.50,
                accrual_fraction_tau: 0.25,
                rfr_forward_volatility_sigma: 0.24 + credit_bump * 2.0,
                discount_factor_p_ti: (-0.75 * r).exp(),
            },
        );
        let ps_coco = ql_fixed_income::evaluate_pennacchi_sundaresan_coco_at1_trigger_valuation(
            &ql_fixed_income::PennacchiSundaresanCocoAt1TriggerSpec {
                current_cet1_ratio_c0: 0.132 - credit_bump,
                trigger_cet1_ratio_c_trig: 0.070,
                cet1_drift_mu_c: 0.008,
                cet1_volatility_sigma_c: 0.19 + credit_bump * 2.0,
                at1_coupon_rate_q: 0.068,
                equity_conversion_share_chi: 0.45,
                write_down_loss_severity_lwd: 0.85,
                risk_free_rate_r: r,
                horizon_years_t: 5.0,
            },
        );
        let pb_cdx = ql_credit::evaluate_pedersen_brigo_cds_index_option_fep_swaption(
            &ql_credit::PedersenBrigoCdsIndexOptionFepSpec {
                current_index_spread_s0: 0.0095 + credit_bump,
                strike_spread_k: 0.0108,
                index_running_coupon_c: 0.0100,
                forward_rpv01_annuity_a0: 4.15,
                spread_volatility_sigma_s: 0.48 + credit_bump * 2.0,
                recovery_rate_r: 0.40,
                risk_free_rate_rf: r,
                option_expiry_years_te: 0.50,
            },
        );
        let bjs_svj = ql_vol::evaluate_broadie_jain_sepp_svj_volatility_swap_convexity(
            &ql_vol::BroadieJainSeppSvjVolatilitySwapSpec {
                initial_variance_v0: vol * vol,
                long_run_variance_theta: 0.036,
                mean_reversion_kappa: 2.10,
                vol_of_vol_xi: 0.45 + credit_bump * 2.0,
                jump_intensity_lambda_j: 0.65 + credit_bump * 5.0,
                jump_mean_mu_j: -0.075,
                jump_std_sigma_j: 0.095,
                maturity_years_t: 1.0,
            },
        );
        let bcp_wwr = ql_xva::evaluate_brigo_capponi_pallavicini_cds_contagion_wwr_cva(
            &ql_xva::BrigoCapponiPallaviciniCdsContagionWwrSpec {
                counterparty_hazard_lambda_c: 0.022 + credit_bump,
                reference_hazard_lambda_r: 0.018 + credit_bump,
                reference_spread_vol_sigma_r: 0.52,
                spread_correlation_rho_rc: 0.55,
                contagion_hazard_jump_eta: 0.65 + credit_bump * 4.0,
                counterparty_lgd: 0.60,
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let rss_store =
            ql_commodities::evaluate_routledge_seppi_spatt_competitive_storage_equilibrium(
                &ql_commodities::RoutledgeSeppiSpattCompetitiveStorageSpec {
                    spot_price_s0: 78.0 + inflation * 40.0,
                    current_inventory_i0: 28.0 - credit_bump * 50.0,
                    critical_stockout_inventory_i_star: 22.0,
                    stockout_convenience_spike_delta: 0.28,
                    storage_cost_rate_c: 0.075,
                    base_spot_volatility_sigma: 0.24 + credit_bump,
                    stockout_vol_amplification_beta: 0.85,
                    maturity_years_t: 1.0,
                },
            );
        let pit_prdc = ql_fx::evaluate_piterbarg_prdc_fx_ir_hybrid_coupon_valuation(
            &ql_fx::PiterbargPrdcHybridSmileCouponSpec {
                spot_fx_s0: 104.0,
                domestic_rate_rd: 0.008,
                foreign_rate_rf: r + 0.012,
                foreign_coupon_gearing_gf: 0.12,
                domestic_coupon_strike_gd: 0.08,
                fx_spot_volatility_sigma_s: 0.115 + credit_bump,
                rate_diff_volatility_sigma_ir: 0.012,
                fx_ir_correlation_rho: 0.35,
                maturity_years_t: 10.0,
            },
        );
        let bgk_bar = ql_exotics::evaluate_broadie_glasserman_kou_discrete_barrier_call(
            &ql_exotics::BroadieGlassermanKouDiscreteBarrierCallSpec {
                spot_price_s0: spot,
                strike_price_k: spot,
                lower_barrier_h: spot * 0.88,
                risk_free_rate_r: r,
                volatility_sigma: vol,
                maturity_years_t: 1.0,
                monitoring_intervals_m: 52.0,
            },
        );
        let bsv_ppp =
            ql_portfolio::evaluate_brandt_santa_clara_valkanov_parametric_portfolio_policy(
                &ql_portfolio::BrandtSantaClaraValkanovPppSpec {
                    value_characteristic_mean_return: 0.045,
                    momentum_characteristic_mean_return: 0.062 - credit_bump,
                    value_characteristic_variance: 0.018,
                    momentum_characteristic_variance: 0.024,
                    value_momentum_covariance: -0.008,
                    benchmark_value_covariance: 0.002,
                    benchmark_momentum_covariance: -0.001,
                    relative_risk_aversion_gamma: 5.0,
                },
            );
        let gss_imp = ql_portfolio::evaluate_gatheral_schied_slynko_transient_impact_execution(
            &ql_portfolio::GatheralSchiedSlynkoTransientImpactSpec {
                order_size_shares_x0: 10.0,
                transient_impact_coeff_eta: 0.085 + credit_bump,
                resilience_decay_rate_rho: 1.60,
                permanent_impact_coeff_gamma: 0.025,
                price_volatility_sigma: vol,
                risk_aversion_lambda: 0.15,
                execution_horizon_hours_t: 4.0,
            },
        );
        let bpq_harq = ql_dml::evaluate_corsi_bollerslev_patton_quaedvlieg_harq_volatility(
            &ql_dml::CorsiBollerslevPattonQuaedvliegHarqSpec {
                intercept_beta0: 0.000018,
                daily_persistence_beta_d: 0.42,
                weekly_persistence_beta_w: 0.34,
                monthly_persistence_beta_m: 0.20,
                quarticity_attenuation_beta_q: -0.16,
                daily_realized_variance_rv_d: 0.00024 + credit_bump * 0.005,
                weekly_realized_variance_rv_w: 0.00018,
                monthly_realized_variance_rv_m: 0.00015,
                realized_quarticity_excess_rq: 0.65 + credit_bump * 5.0,
            },
        );
        let ll_cay = ql_models::evaluate_lettau_ludvigson_cay_conditional_equity_premium(
            &ql_models::LettauLudvigsonCayConditionalRiskPremiumSpec {
                log_consumption_c: 10.24 + inflation * 0.5,
                log_asset_wealth_a: 10.05 - credit_bump * 2.0,
                log_labor_income_y: 10.28,
                asset_wealth_cointegration_share_omega: 0.32,
                baseline_risk_aversion_gamma0: 12.0,
                cay_conditioning_sensitivity_gamma1: 185.0,
                consumption_volatility_sigma_c: 0.018,
                equity_return_volatility_sigma_m: vol,
            },
        );

        lm_rfr.in_arrears_rfr_caplet_pv * 10.0
            + ps_coco.coco_at1_fair_price_pct * 0.1
            + pb_cdx.cds_index_payer_swaption_pv_pct * 10.0
            + bjs_svj.svj_fair_volatility_strike * 10.0
            + bcp_wwr.total_contagion_wwr_cva_pct * 10.0
            + rss_store.equilibrium_futures_price * 0.1
            + pit_prdc.cancelable_prdc_net_value
            + bgk_bar.discrete_bgk_down_and_out_call
            + bsv_ppp.certainty_equivalent_utility_gain_bps * 0.01
            + gss_imp.expected_implementation_shortfall_cost
            + bpq_harq.annualized_harq_volatility_forecast * 10.0
            + ll_cay.conditional_equity_risk_premium * 100.0
    }
    pub fn lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let cp_cycle = ql_curves::evaluate_cieslak_povala_trend_inflation_bond_cycle(
            &ql_curves::CieslakPovalaTrendInflationCycleSpec {
                short_rate_y1: r,
                intermediate_yield_y5: r + 0.004,
                long_yield_y10: r + 0.008,
                discounted_cpi_trend_tau: inflation,
                short_rate_inflation_loading_b1: 1.05,
                long_yield_inflation_loading_b10: 0.92,
                real_rate_cycle_return_loading_gamma_c: 0.28,
                short_rate_cycle_return_loading_gamma_1: 0.14,
            },
        );
        let gkv_mbs = ql_fixed_income::evaluate_gabaix_krishnamurthy_vigneron_mbs_intermediary_oas(
            &ql_fixed_income::GabaixKrishnamurthyVigneronMbsOasSpec {
                mortgage_coupon_wac: 0.055,
                current_refi_rate_r: r + 0.015,
                prepayment_s_curve_slope_beta: 25.0,
                rate_volatility_sigma_r: 0.012 + credit_bump * 0.2,
                aggregate_mbs_outstanding_share_w: 0.65,
                intermediary_risk_aversion_gamma: 18.0 + credit_bump * 50.0,
                pass_through_duration_d0: 5.4,
                liquidity_hedging_cost_lambda: 3.5,
            },
        );
        let egg_hawkes = ql_credit::evaluate_errais_giesecke_goldberg_hawkes_contagion_credit_loss(
            &ql_credit::ErraisGieseckeGoldbergHawkesCreditSpec {
                initial_default_intensity_lambda0: 0.035 + credit_bump,
                long_run_baseline_intensity_mu: 0.025 + credit_bump * 0.5,
                mean_reversion_speed_kappa: 1.40,
                self_exciting_contagion_impact_delta: 0.95,
                expected_loss_per_default_ell: 0.60,
                horizon_years_t: 5.0,
                tranche_detachment_fraction_d: 0.12,
            },
        );
        let mb_svix = ql_vol::evaluate_martin_bondarenko_simple_variance_svix_lower_bound(
            &ql_vol::MartinBondarenkoSimpleVarianceSvixSpec {
                forward_price_f0: spot * 45.0,
                risk_free_rate_r: r,
                atm_volatility_sigma: vol,
                downside_skew_slope_psi: -0.18 - credit_bump * 2.0,
                crash_jump_intensity_lambda: 0.45 + credit_bump * 5.0,
                crash_log_jump_mean_mu_j: -0.14,
                relative_risk_aversion_gamma: 2.5,
                horizon_years_t: 0.25,
            },
        );
        let cbc_bsde = ql_xva::evaluate_crepey_bichuch_capponi_nonlinear_xva_bsde(
            &ql_xva::CrepeyBichuchCapponiNonlinearXvaBsdeSpec {
                clean_derivative_mtm_v0: 12.5,
                collateralization_fraction_alpha: 0.35,
                unsecured_borrowing_spread_s_b: 0.0145 + credit_bump,
                unsecured_lending_spread_s_l: 0.0040,
                counterparty_default_spread_s_c: 0.0120 + credit_bump,
                hedger_default_spread_s_h: 0.0065,
                exposure_volatility_sigma_v: 8.5,
                maturity_years_t: 5.0,
            },
        );
        let ts_usv =
            ql_commodities::evaluate_trolle_schwartz_unspanned_stochastic_volatility_commodity(
                &ql_commodities::TrolleSchwartzUnspannedStochasticVolatilitySpec {
                    futures_price_f0: 76.0 + inflation * 50.0,
                    forward_cost_of_carry_y0: -0.035,
                    carry_decay_kappa_y: 1.25,
                    initial_stochastic_variance_v0: 0.095 + credit_bump * 2.0,
                    long_run_variance_theta_v: 0.080,
                    variance_mean_reversion_kappa_v: 1.80,
                    vol_of_vol_sigma_v: 0.42,
                    option_expiry_years_t: 0.50,
                },
            );
        let iss_cip = ql_fx::evaluate_ivashina_scharfstein_stein_dollar_funding_cip_basis(
            &ql_fx::IvashinaScharfsteinSteinDollarFundingCipSpec {
                direct_usd_commercial_paper_rate_r_usd: r + 0.015,
                foreign_domestic_deposit_rate_r_for: r - 0.005,
                wholesale_usd_deposit_outflow_shock_delta: 0.25 + credit_bump * 4.0,
                dollar_loan_book_rigidity_theta: 0.55,
                fx_swap_dealer_balance_sheet_cost_eta: 0.018,
                spot_fx_s0: 1.08,
                credit_worthiness_haircut_alpha: 0.012,
                tenor_years_t: 0.50,
            },
        );
        let ww_cliq = ql_exotics::evaluate_wilmott_windcliff_cliquet_forward_skew_valuation(
            &ql_exotics::WilmottWindcliffCliquetForwardSkewSpec {
                notional_amount: 100.0,
                local_floor_f_loc: 0.0,
                local_cap_c_loc: 0.08,
                global_floor_f_glob: 0.16,
                risk_free_rate_r: r,
                atm_reset_volatility_sigma: vol,
                upper_cap_skew_slope_beta: -0.35 - credit_bump * 2.0,
                reset_periods_m: 5.0,
            },
        );
        let lw_nl = ql_portfolio::evaluate_ledoit_wolf_analytical_nonlinear_eigenvalue_shrinkage(
            &ql_portfolio::LedoitWolfNonlinearEigenvalueShrinkageSpec {
                sample_eigenvalue_lambda: 0.42 + credit_bump * 2.0,
                concentration_ratio_c: 0.35,
                bulk_spectral_mean_bar_lambda: 1.00,
                local_spectral_density_f: 0.68,
                hilbert_transform_stieltjes_h: 0.45,
                tail_regularization_weight_xi: 0.08,
            },
        );
        let ab_covar = ql_portfolio::evaluate_adrian_brunnermeier_systemic_delta_covar(
            &ql_portfolio::AdrianBrunnermeierSystemicCoVarSpec {
                institution_mean_loss_mu_i: 0.004,
                institution_volatility_sigma_i: 0.032 + credit_bump,
                system_mean_loss_mu_s: 0.003,
                system_volatility_sigma_s: 0.024 + credit_bump,
                normal_correlation_rho: 0.58,
                distress_quantile_slope_delta_beta: 0.22,
                leverage_maturity_mismatch_m: 1.45 + credit_bump * 10.0,
                tail_quantile_z_q: 2.3263,
            },
        );
        let btz_vrp = ql_dml::evaluate_bollerslev_tauchen_zhou_variance_risk_premium_predictability(
            &ql_dml::BollerslevTauchenZhouVarianceRiskPremiumSpec {
                model_free_implied_variance_iv: vol * vol + credit_bump * 0.5,
                high_frequency_realized_variance_rv: 0.024,
                variance_persistence_rho_v: 0.62,
                vol_of_vol_persistence_rho_q: 0.82,
                vol_of_vol_shock_variance_varphi: 0.18,
                consumption_volatility_sensitivity_a_q: 2.40,
                recursive_risk_aversion_gamma: 6.5,
                forecast_horizon_quarters_h: 2.0,
            },
        );
        let cv_beta = ql_models::evaluate_campbell_vuolteenaho_bad_beta_good_beta_icapm(
            &ql_models::CampbellVuolteenahoBadGoodBetaIcapmSpec {
                cash_flow_news_covariance_cov_cf: 0.0115 + credit_bump * 0.2,
                discount_rate_news_covariance_cov_dr: 0.0145,
                market_return_variance_var_m: vol * vol,
                investor_relative_risk_aversion_gamma: 7.5,
                log_linearization_discount_rho: 0.95,
                return_forecasting_var_persistence_phi: 0.88,
                idiosyncratic_alpha_adjustment_xi: 0.015,
            },
        );

        cp_cycle.expected_10y_bond_excess_return * 100.0
            + gkv_mbs.equilibrium_mbs_total_oas_bps * 0.1
            + egg_hawkes.contagion_adjusted_tranche_risk_charge * 10.0
            + mb_svix.martin_generalized_equity_premium_bound * 10.0
            + cbc_bsde.nonlinear_bsde_total_xva
            + ts_usv.usv_atm_futures_option_value
            + iss_cip.segmented_fx_forward_rate * 10.0
            + ww_cliq.cliquet_forward_skew_fair_value * 0.1
            + lw_nl.nonlinear_shrunk_eigenvalue * 10.0
            + ab_covar.system_covar_in_distress * 100.0
            + btz_vrp.predicted_annualized_excess_return * 10.0
            + cv_beta.icapm_two_beta_required_risk_premium * 100.0
    }
    pub fn cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let bm_fra = ql_curves::evaluate_bianchetti_mercurio_multi_curve_fra_quanto_convexity(
            &ql_curves::BianchettiMercurioMultiCurveFraQuantoSpec {
                ois_forward_rate_f_d: r,
                ibor_funding_spread_s0: 0.0045 + credit_bump,
                fixed_fra_rate_k: 0.040,
                ois_forward_volatility_sigma_d: 0.22,
                spread_volatility_sigma_s: 0.35 + credit_bump * 2.0,
                ois_spread_correlation_rho: -0.25,
                fixing_time_years_t1: 1.5,
                accrual_fraction_tau: 0.50,
            },
        );
        let aab_liq = ql_fixed_income::evaluate_acharya_amihud_bharath_liquidity_stress_bond_spread(
            &ql_fixed_income::AcharyaAmihudBharathLiquidityStressBondSpec {
                structural_default_spread_s_def: 0.018 + credit_bump,
                bond_illiquidity_level_illiq: 0.012 + credit_bump * 0.5,
                normal_regime_liquidity_beta_1: 0.25,
                stress_regime_liquidity_beta_2: 0.85,
                stress_regime_probability_pi_s: 0.30 + credit_bump * 5.0,
                market_illiquidity_risk_price_lambda_l: 0.015,
                flight_to_quality_treasury_wedge_psi: 0.004,
                modified_duration_d: 6.2,
            },
        );
        let ba_cds = ql_credit::evaluate_brigo_alfonsi_correlated_cir_stochastic_intensity_cds(
            &ql_credit::BrigoAlfonsiCorrelatedCirCdsSpec {
                initial_short_rate_r0: r,
                initial_default_intensity_lambda0: 0.022 + credit_bump,
                intensity_mean_reversion_kappa: 0.65,
                intensity_long_run_theta: 0.028 + credit_bump * 0.5,
                intensity_volatility_sigma_lambda: 0.18,
                rate_volatility_sigma_r: 0.08,
                rate_intensity_correlation_rho: -0.35,
                loss_given_default_lgd: 0.60,
            },
        );
        let fuk_skew = ql_vol::evaluate_fukasawa_asymptotic_moneyness_skew_varswap(
            &ql_vol::FukasawaAsymptoticImpliedSkewVarswapSpec {
                atm_implied_volatility_sigma0: vol,
                standardized_moneyness_skew_a1: -0.035 - credit_bump,
                standardized_moneyness_curvature_a2: 0.018,
                third_order_wing_asymmetry_a3: -0.006,
                risk_free_rate_r: r,
                maturity_years_t: 0.50,
                target_standardized_moneyness_z: -1.25,
            },
        );
        let ads_fva = ql_xva::evaluate_andersen_duffie_song_debt_overhang_fva(
            &ql_xva::AndersenDuffieSongShareholderFvaOverhangSpec {
                uncollateralized_funding_need_phi: 100.0,
                dealer_credit_spread_s_d: 0.0125 + credit_bump,
                upfront_marginal_debt_issuance_share_q: 0.80,
                existing_creditor_recovery_dilution_alpha: 0.25,
                asset_swap_covariance_beta_v: 0.08,
                risk_free_rate_r: r,
                marginal_capital_return_hurdle_eta: 0.018,
                maturity_years_t: 5.0,
            },
        );
        let pir_sprd = ql_commodities::evaluate_pirrong_structural_inventory_skew_spread_option(
            &ql_commodities::PirrongStructuralInventorySkewSpreadSpec {
                prompt_futures_f1: 82.0 + inflation * 50.0,
                deferred_futures_f2: 79.5,
                spread_strike_k: 2.0,
                normalized_inventory_scarcity_z_inv: 1.15 + credit_bump * 5.0,
                baseline_prompt_vol_sigma1: 0.28,
                deferred_vol_sigma2: 0.22,
                baseline_correlation_rho0: 0.92,
                scarcity_decorrelation_sensitivity_psi: 0.18,
            },
        );
        let lc_dnt = ql_fx::evaluate_lipton_castagna_fx_double_no_touch_vanna_volga(
            &ql_fx::LiptonCastagnaFxDoubleNoTouchVannaVolgaSpec {
                spot_fx_s0: 1.092,
                lower_barrier_l: 1.030,
                upper_barrier_u: 1.150,
                domestic_rate_rd: r + 0.010,
                foreign_rate_rf: r - 0.005,
                atm_volatility_sigma: 0.088 + credit_bump,
                risk_reversal_25d_rr: -0.012,
                butterfly_25d_bf: 0.0045 + credit_bump * 0.2,
            },
        );
        let lin_step = ql_exotics::evaluate_linetsky_occupation_time_step_barrier_option(
            &ql_exotics::LinetskyStepBarrierOccupationSpec {
                spot_price_s0: spot,
                strike_price_k: spot,
                lower_step_barrier_h: spot * 0.88,
                risk_free_rate_r: r,
                volatility_sigma: vol,
                occupation_killing_rate_rho: 8.0 + credit_bump * 20.0,
                maturity_years_t: 1.0,
            },
        );
        let ky_dem = ql_portfolio::evaluate_koijen_yogo_characteristic_demand_system_elasticity(
            &ql_portfolio::KoijenYogoDemandSystemAssetPricingSpec {
                log_market_to_book_me_be: 0.65,
                beta_loading_beta_m: 1.05,
                profitability_characteristic_prof: 0.32 - credit_bump,
                institutional_price_elasticity_coeff_beta0: -0.45,
                characteristic_tilt_strength_beta_x: 1.20,
                outside_asset_portfolio_weight_w0: 0.22,
                active_institutional_aum_share_theta: 0.68,
                inelastic_passive_flow_shock_f: 0.025 - credit_bump,
            },
        );
        let ko_inv = ql_portfolio::evaluate_kyle_obizhaeva_microstructure_invariance_execution(
            &ql_portfolio::KyleObizhaevaMicrostructureInvarianceSpec {
                order_size_shares_q: 25_000.0,
                daily_trading_volume_v: 1_200_000.0,
                stock_price_p: spot,
                daily_return_volatility_sigma: 0.018 + credit_bump * 0.5,
                invariance_dollar_bet_var_i_star: 25_000.0,
                linear_impact_calibration_c_lambda: 0.32,
                bid_ask_spread_calibration_c_s: 0.08,
            },
        );
        let gkx_r2 = ql_dml::evaluate_gu_kelly_xiu_nonlinear_characteristic_oos_r2(
            &ql_dml::GuKellyXiuPredictiveOosR2InteractionSpec {
                momentum_signal_z_mom: 0.75,
                value_signal_z_val: 0.55,
                macro_equity_premium_spread_x_m: 0.028 + credit_bump,
                linear_feature_weight_w1: 0.42,
                nonlinear_interaction_weight_w12: 0.35,
                l2_ridge_penalization_lambda: 0.15,
                realized_excess_return_r_next: 0.019,
                cross_sectional_return_dispersion_s2: 0.0085,
            },
        );
        let dss_dep = ql_models::evaluate_drechsler_savov_schnabl_deposit_channel_liquidity_premium(
            &ql_models::DrechslerSavovSchnablDepositChannelSpec {
                fed_funds_policy_rate_f: r + 0.015,
                bank_branch_herfindahl_hhi: 0.24,
                household_cash_substitution_elasticity_eps: 0.55,
                bank_liquidity_management_cost_c: 0.0045,
                wholesale_funding_sensitivity_delta: 0.35,
                loan_illiquidity_spread_alpha: 1.80,
                aggregate_liquidity_demand_kappa: 1.25,
                policy_rate_hike_shock_df: 0.010 + credit_bump,
            },
        );

        bm_fra.multi_curve_payer_fra_pv * 1000.0
            + aab_liq.stress_adjusted_corporate_bond_spread * 100.0
            + ba_cds.correlated_ssrd_par_cds_spread_bps * 0.01
            + fuk_skew.fukasawa_variance_swap_fair_strike * 100.0
            + ads_fva.shareholder_indifference_fva
            + pir_sprd.pirrong_calendar_spread_call_value
            + lc_dnt.vanna_volga_dnt_fair_value * 10.0
            + lin_step.step_barrier_down_and_out_call
            + ky_dem.equilibrium_log_price_impact * 100.0
            + ko_inv.invariance_total_dollar_execution_cost * 1e-3
            + gkx_r2.out_of_sample_predictive_r2 * 10.0
            + dss_dep.equilibrium_liquidity_premium_wedge * 100.0
    }
    pub fn bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let hs_dur = ql_curves::evaluate_hanson_stein_greenwood_vayanos_duration_term_premium(
            &ql_curves::HansonSteinGreenwoodVayanosDurationSupplySpec {
                short_real_rate_r1: r * 0.5,
                expected_future_short_rate_er2: r * 0.55,
                net_duration_supply_s: 4.5,
                mbs_convexity_supply_shock_m: 0.85 + credit_bump * 10.0,
                arbitrageur_risk_aversion_gamma: 22.0,
                short_rate_volatility_sigma_r: 0.014 + credit_bump * 0.2,
                yield_oriented_investor_share_alpha: 0.45,
                reaching_for_yield_sensitivity_theta: 18.0,
            },
        );
        let fll_tips = ql_fixed_income::evaluate_fleckenstein_longstaff_lustig_tips_treasury_puzzle(
            &ql_fixed_income::FleckensteinLongstaffLustigTipsTreasuryPuzzleSpec {
                nominal_treasury_yield_y_nom: r + 0.010,
                tips_real_yield_y_tips: r * 0.65 + credit_bump,
                zero_coupon_inflation_swap_rate_s_cpi: inflation,
                tips_deflation_floor_strike_k_floor: 0.0,
                cpi_inflation_volatility_sigma_pi: 0.015 + credit_bump * 0.5,
                repo_specialness_differential_delta_repo: 0.0012,
                dealer_balance_sheet_leverage_cost_eta: 0.35,
                maturity_years_t: 5.0,
            },
        );
        let cdgy_cdo = ql_credit::evaluate_collin_dufresne_goldstein_yang_super_senior_cdo(
            &ql_credit::CollinDufresneGoldsteinYangSuperSeniorCdoSpec {
                idiosyncratic_hazard_rate_lambda_i: 0.015 + credit_bump,
                catastrophic_macro_jump_intensity_lambda_m: 0.018 + credit_bump,
                macro_crisis_conditional_default_prob_p_m: 0.36 + credit_bump * 2.0,
                loss_given_default_lgd: 0.60,
                super_senior_attachment_a: 0.15,
                super_senior_detachment_d: 0.30,
                risk_free_rate_r: r,
                maturity_years_t: 5.0,
            },
        );
        let tod_vol = ql_vol::evaluate_todorov_characteristic_option_spot_volatility(
            &ql_vol::TodorovOptionImpliedSpotVolatilityJumpSpec {
                short_tenor_years_tau: 0.020,
                characteristic_wavenumber_u: 2.5,
                empirical_log_characteristic_modulus_l_u: 0.125 + vol * 0.1 + credit_bump,
                high_frequency_wavenumber_ratio_k: 1.6,
                jump_activity_blumenthal_getoor_beta: 1.25,
                tail_jump_intensity_scale_c_j: 0.035 + credit_bump,
                risk_free_rate_r: r,
            },
        );
        let acc_kva = ql_xva::evaluate_albanese_caenazzo_crepey_bilateral_kva_mva(
            &ql_xva::AlbaneseCaenazzoCrepeyBilateralKvaMvaSpec {
                regulatory_capital_ead: 120.0 + credit_bump * 200.0,
                counterparty_risk_weight_rw: 0.08 + credit_bump,
                cva_capital_charge_k_cva: 3.5,
                shareholder_hurdle_rate_h: 0.105,
                initial_margin_posted_im0: 45.0,
                im_funding_spread_s_im: 0.0085 + credit_bump,
                capital_decay_rate_lambda_k: 0.18,
                horizon_years_t: 7.0,
            },
        );
        let alr_hp =
            ql_commodities::evaluate_acharya_lochstoer_ramadorai_commodity_hedging_pressure(
                &ql_commodities::AcharyaLochstoerRamadoraiHedgingPressureSpec {
                    expected_future_spot_e_s1: 80.0 + inflation * 50.0,
                    spot_price_volatility_sigma_s: 14.0,
                    producer_default_prob_zm_score: 0.45 + credit_bump * 5.0,
                    producer_inventory_endowment_q_p: 2.5,
                    producer_distress_aversion_gamma_p: 1.8,
                    speculator_capital_risk_aversion_gamma_s: 0.65,
                    speculator_wealth_constraint_tightness_w: 0.55 + credit_bump * 4.0,
                    risk_free_rate_r: r,
                },
            );
        let bnp_fx = ql_fx::evaluate_brunnermeier_nagel_pedersen_fx_carry_crash_skew(
            &ql_fx::BrunnermeierNagelPedersenFxCarryCrashSpec {
                interest_rate_differential_idiff: 0.045,
                speculator_carry_position_size_x: 1.40,
                funding_liquidity_ted_spread_s_ted: 0.018 + credit_bump,
                margin_spiral_amplification_phi: 4.5,
                baseline_fx_volatility_sigma_fx: 0.115 + credit_bump,
                crash_arrival_sensitivity_kappa: 2.2,
                risk_reversal_put_skew_cost_c_rr: 0.35,
                horizon_years_t: 0.50,
            },
        );
        let by_asian = ql_exotics::evaluate_broadie_yamamoto_fast_gauss_discrete_asian_option(
            &ql_exotics::BroadieYamamotoFastGaussDiscreteAsianSpec {
                spot_price_s0: spot,
                strike_price_k: spot,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                discrete_fixings_m: 12.0,
                maturity_years_t: 1.0,
            },
        );
        let kps_ipca = ql_portfolio::evaluate_kelly_pruitt_su_instrumented_pca_dynamic_beta(
            &ql_portfolio::KellyPruittSuInstrumentedPcaDynamicBetaSpec {
                characteristic_value_z1: 0.65,
                characteristic_momentum_z2: 0.80 - credit_bump,
                characteristic_investment_z3: 0.40,
                gamma_beta_loading_1: 0.55,
                gamma_beta_loading_2: 0.72,
                gamma_alpha_anomaly_intercept: 0.015,
                latent_factor_realized_return_f: 0.068 - credit_bump,
                cross_sectional_idiosyncratic_var_sigma2: 0.042,
            },
        );
        let be_srisk = ql_portfolio::evaluate_brownlees_engle_srisk_capital_shortfall(
            &ql_portfolio::BrownleesEngleSystemicSriskCapitalShortfallSpec {
                book_value_of_debt_d: 950.0,
                market_value_of_equity_w: spot * 0.85,
                prudential_capital_ratio_k: 0.08,
                dynamic_conditional_beta_rho: 1.35 + credit_bump * 5.0,
                systemic_crisis_market_decline_c: 0.40 + credit_bump * 2.0,
                idiosyncratic_tail_amplification_theta: 0.18,
                funding_run_haircut_psi: 0.05,
            },
        );
        let abdl_rv = ql_dml::evaluate_andersen_bollerslev_diebold_labys_fractional_rv(
            &ql_dml::AndersenBollerslevDieboldLabysArfimaRvSpec {
                current_realized_variance_rv0: 0.00022 + credit_bump * 0.002,
                lag1_realized_variance_rv1: 0.00018,
                long_run_mean_log_rv_mu: -8.85,
                fractional_differencing_d: 0.40,
                short_run_ar1_phi: 0.18,
                log_rv_innovation_std_sigma_u: 0.32,
                forecast_horizon_days_h: 5.0,
            },
        );
        let fg_uip = ql_models::evaluate_farhi_gabaix_rare_disaster_exchange_rate_uip(
            &ql_models::FarhiGabaixRareDisasterExchangeRatePuzzleSpec {
                global_disaster_probability_p: 0.036 + credit_bump,
                world_disaster_recovery_b_star: 2.40,
                home_country_resilience_h: 0.15 - credit_bump,
                foreign_country_resilience_h_star: 0.45,
                resilience_mean_reversion_phi_h: 0.18,
                domestic_productivity_growth_g: 0.018,
                real_discount_rate_r_bar: r,
            },
        );

        hs_dur.equilibrium_long_real_yield * 100.0
            + fll_tips.tips_treasury_arbitrage_wedge_per_100
            + cdgy_cdo.super_senior_tranche_par_spread_bps * 0.01
            + tod_vol.jump_robust_spot_volatility * 10.0
            + acc_kva.total_bilateral_kva_plus_mva
            + alr_hp.equilibrium_futures_price * 0.1
            + bnp_fx.crash_adjusted_expected_carry_return * 100.0
            + by_asian.fast_gauss_discrete_asian_call
            + kps_ipca.ipca_precision_weighted_expected_return * 100.0
            + be_srisk.conditional_capital_shortfall_srisk * 0.1
            + abdl_rv.annualized_arfima_volatility_forecast * 10.0
            + fg_uip.expected_disaster_fx_carry_excess_return * 100.0
    }
    pub fn hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }
    pub fn klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.05;
        let spot = m.eq_spot.max(50.0);
        let vol = m.eq_vol.max(0.15);

        let ks_sofr = ql_curves::evaluate_klingler_syrstad_sofr_repo_balance_sheet_spike(
            &ql_curves::KlinglerSyrstadSofrRepoSpecialnessSpikeSpec {
                ois_policy_rate_r_ois: r + 0.015,
                treasury_collateral_specialness_s_spec: 0.0018,
                dealer_treasury_inventory_congestion_q: 0.65 + credit_bump * 5.0,
                quarter_end_balance_sheet_spike_a_qe: 0.0085 + credit_bump * 0.5,
                reporting_date_weight_omega_qe: 0.08,
                bank_credit_sensitive_spread_s_cr: 0.0032 + credit_bump,
                repo_rate_volatility_sigma_repo: 0.012,
                tenor_years_tau: 0.25,
            },
        );
        let lmn_basis = ql_fixed_income::evaluate_longstaff_mithal_neis_cds_bond_liquidity_basis(
            &ql_fixed_income::LongstaffMithalNeisCdsBondLiquidityBasisSpec {
                corporate_bond_yield_y_corp: r + 0.025 + credit_bump * 1.5,
                treasury_benchmark_yield_y_tsy: r,
                cds_default_spread_s_cds: 0.0115 + credit_bump,
                corporate_coupon_rate_c: 0.055,
                state_tax_differential_tau: 0.065,
                bond_bid_ask_spread_s_ba: 0.0048 + credit_bump * 0.5,
                issue_turnover_liquidity_factor_l: 1.25,
                modified_duration_d: 6.0,
            },
        );
        let ps_sov = ql_credit::evaluate_pan_singleton_sovereign_cds_distress_risk_premium(
            &ql_credit::PanSingletonSovereignCdsDistressRiskSpec {
                physical_hazard_rate_lambda0_p: 0.018 + credit_bump,
                physical_mean_reversion_kappa_p: 0.55,
                log_intensity_volatility_sigma: 0.42 + credit_bump * 2.0,
                distress_risk_price_eta: 0.45,
                risk_neutral_loss_given_default_lgd: 0.75,
                short_cds_tenor_t1: 1.0,
                long_cds_tenor_t2: 5.0,
                risk_free_rate_r: r,
            },
        );
        let ms_smile = ql_vol::evaluate_medvedev_scaillet_short_maturity_smile_expansion(
            &ql_vol::MedvedevScailletShortMaturitySmileExpansionSpec {
                spot_volatility_sigma0: vol,
                vol_of_vol_nu: 0.55 + credit_bump * 2.0,
                spot_vol_correlation_rho: -0.68,
                volatility_drift_mu_v: 0.015,
                normalized_moneyness_theta: -1.20,
                maturity_years_tau: 0.25,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
            },
        );
        let dz_ccp = ql_xva::evaluate_duffie_zhu_central_clearing_netting_efficiency(
            &ql_xva::DuffieZhuCentralClearingNettingEfficiencySpec {
                clearing_member_count_n: 12.0,
                asset_class_count_k: 4.0,
                cleared_class_share_alpha: 0.85,
                bilateral_exposure_volatility_sigma: 25.0 + credit_bump * 100.0,
                cross_asset_correlation_rho: 0.20,
                default_fund_mutualization_weight_phi: 0.35,
            },
        );
        let ghr_inv = ql_commodities::evaluate_gorton_hayashi_rouwenhorst_commodity_inventory_basis(
            &ql_commodities::GortonHayashiRouwenhorstInventoryBasisSpec {
                normalized_inventory_level_i: 0.72 - inflation * 2.0,
                convenience_yield_scarcity_alpha: 0.085,
                stockout_convexity_exponent_beta: 1.65,
                interest_plus_storage_cost_rc: r + 0.04,
                prior_futures_momentum_signal_m: 0.14,
                inventory_risk_premium_loading_lambda: 0.80,
                spot_volatility_sigma_s: 0.26 + credit_bump,
            },
        );
        let drs_vrp = ql_fx::evaluate_della_corte_ramadorai_sarno_fx_volatility_premium(
            &ql_fx::DellaCorteRamadoraiSarnoFxVolatilityPremiumSpec {
                realized_fx_volatility_rv: 0.118 + credit_bump,
                model_free_implied_fx_vol_iv: 0.092 + credit_bump * 0.5,
                butterfly_25d_smile_convexity_bf: 0.0038,
                forward_discount_carry_idiff: 0.025,
                vrp_spot_predictability_beta: 0.65,
                volatility_insurance_hedging_gamma: 0.40,
                horizon_years_tau: 0.50,
            },
        );
        let pel_db = ql_exotics::evaluate_pelsser_laplace_double_barrier_knockout_option(
            &ql_exotics::PelsserLaplaceDoubleBarrierKnockoutSpec {
                spot_price_s0: spot,
                lower_barrier_l: spot * 0.82,
                upper_barrier_u: spot * 1.22,
                knockout_rebate_r_reb: 15.0,
                risk_free_rate_r: r,
                dividend_yield_q: 0.015,
                volatility_sigma: vol,
                maturity_years_t: 0.75,
            },
        );
        let mm_vmp = ql_portfolio::evaluate_moreira_muir_volatility_managed_portfolio_alpha(
            &ql_portfolio::MoreiraMuirVolatilityManagedPortfolioSpec {
                unconditional_factor_expected_return_mu: 0.075,
                unconditional_factor_variance_bar_var: 0.032,
                lagged_realized_variance_rv_prev: 0.022 + credit_bump * 0.5,
                volatility_scaling_constant_c: 0.030,
                return_variance_predictability_slope_phi: -0.05,
                turnover_transaction_cost_kappa: 0.0020,
                variance_persistence_rho: 0.72,
            },
        );
        let bcs_fba = ql_portfolio::evaluate_budish_cramton_shim_frequent_batch_auction_sniping(
            &ql_portfolio::BudishCramtonShimFrequentBatchAuctionSpec {
                fundamental_jump_intensity_lambda_j: 1.8 + credit_bump * 10.0,
                investor_arrival_intensity_lambda_i: 8.5,
                mean_jump_size_sigma_j: 0.045,
                competing_hft_count_n: 6.0,
                speed_technology_cost_c_speed: 0.012,
                reaction_latency_delta: 0.001,
                batch_auction_interval_tau: 0.100,
            },
        );
        let amz_tsrv = ql_dml::evaluate_ait_sahalia_mykland_zhang_two_scales_rv(
            &ql_dml::AitSahaliaMyklandZhangTwoScalesRvSpec {
                full_grid_realized_variance_rv_all: 0.00052 + credit_bump * 0.005,
                subsampled_average_rv_avg: 0.00021 + credit_bump * 0.002,
                total_intraday_ticks_n: 4800.0,
                subsampling_grids_k: 30.0,
                integrated_quarticity_iq: 1.5e-7,
            },
        );
        let gk_qe = ql_models::evaluate_gertler_karadi_intermediary_credit_policy(
            &ql_models::GertlerKaradiUnconventionalMonetaryPolicySpec {
                risk_free_deposit_rate_r: r,
                moral_hazard_divertible_fraction_lambda: 0.38,
                franchise_net_worth_multiplier_nu: 1.55 - credit_bump * 5.0,
                excess_return_shadow_value_mu: 0.042,
                central_bank_intermediation_share_psi: 0.18 + credit_bump * 2.0,
                central_bank_efficiency_cost_tau: 0.0015,
                capital_quality_shock_xi: 0.25 + credit_bump * 5.0,
            },
        );

        ks_sofr.credit_vs_sofr_term_basis_bps * 0.1
            + lmn_basis.adjusted_cds_bond_basis_bps * 0.1
            + ps_sov.long_tenor_sovereign_cds_bps * 0.05
            + ms_smile.asymptotic_implied_volatility * 10.0
            + dz_ccp.unified_multilateral_ccp_epe * 0.1
            + ghr_inv.expected_commodity_futures_excess_return * 10.0
            + drs_vrp.vrp_strategy_expected_excess_return * 100.0
            + pel_db.pelsser_double_barrier_total_value * 0.1
            + mm_vmp.volatility_managed_spanning_alpha * 10.0
            + bcs_fba.sniping_latency_tax_elimination_gain * 100.0
            + amz_tsrv.annualized_tsrv_volatility * 10.0
            + gk_qe.equilibrium_credit_spread_bps * 0.05
    }
    pub fn klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric(
        &self,
    ) -> f64 {
        let bump = match self.scenario {
            ScenarioKind::CreditWideningPlus100Bp => 0.01,
            ScenarioKind::StagflationPlus150Bp => 0.005,
            _ => 0.0,
        };
        Self::klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let inflation = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump + inflation * 0.01;
        let f = m.eq_spot.max(1.0);
        let vols = m.eq_vol.max(0.1);
        let c = ql_curves::evaluate_brigo_mercurio_rebonato_swaption_and_cms(
            &ql_curves::BrigoMercurioRebonatoCmsSpec {
                fixing_tenors: vec![1.0, 2.0],
                accrual_taus: vec![1.0, 1.0],
                forward_rates: vec![r, r + 0.005],
                discount_factors: vec![(-r).exp(), (-2.0 * r).exp()],
                forward_vols: vec![vols, vols],
                corr_decay_beta: 0.1,
                option_expiry: 1.0,
                swaption_strike: r,
            },
        );
        let cr = ql_credit::evaluate_duffie_singleton_leland_toft_capital_structure(
            &ql_credit::DuffieSingletonLelandToftSpec {
                unlevered_asset_value: f,
                annual_coupon_c: 0.04 * f,
                debt_principal_p: 0.5 * f,
                rollover_maturity_t: 5.0,
                corporate_tax_rate: 0.25,
                bankruptcy_cost_alpha: 0.2,
                asset_payout_rate: 0.02,
                risk_free_rate: r.max(0.001),
                asset_volatility: vols,
            },
        );
        let ca = ql_fixed_income::evaluate_calamos_despiegeleer_ascot_and_gamma_scalp(
            &ql_fixed_income::CalamosDeSpiegeleerAscotGammaSpec {
                spot_price: f,
                reference_spot: f,
                bond_par: 100.0,
                coupon_rate: 0.04,
                conversion_ratio: 1.0,
                implied_vol: vols,
                realized_vol: vols * 1.1,
                risk_free_rate: r.max(0.001),
                dividend_yield: 0.01,
                base_hazard_rate: 0.02,
                equity_credit_power_p: 1.2,
                recovery_rate: 0.4,
                maturity_years: 3.0,
                holding_period_years: 1.0,
            },
        );
        let sp = ql_fixed_income::evaluate_tuckman_specified_pool_payup_and_dollar_roll(
            &ql_fixed_income::TuckmanSpecifiedPoolDollarRollSpec {
                net_coupon_rate: 0.03,
                discount_yield: r.max(0.001),
                effective_maturity_years: 5.0,
                ctd_tba_cpr: 0.08,
                llb_pool_cpr_multiplier: 0.7,
                tba_front_month_price: 100.0,
                tba_back_month_price: 99.8,
                gc_repo_rate: r,
            },
        );
        let og = ql_portfolio::evaluate_ohara_glosten_milgrom_sequential_spread(
            &ql_portfolio::OharaGlostenMilgromSequentialSpec {
                low_value_vl: f * 0.9,
                high_value_vh: f * 1.1,
                initial_prior_p0: 0.5,
                informed_trader_fraction_mu: 0.3,
                informed_signal_accuracy_theta: 0.8,
                dealer_inventory_gamma: 0.01,
                signed_order_flow: vec![0.2, -0.1, 0.3],
            },
        );
        let rie = ql_portfolio::evaluate_potters_bouchaud_stieltjes_rie_shrinkage(
            &ql_portfolio::PottersBouchaudStieltjesRieSpec {
                sample_eigenvalues: vec![0.5 + credit_bump, 1.0 + r * 0.01, 2.0 + vols * 0.01],
                aspect_ratio_q: (0.5 + credit_bump).min(0.9),
                stieltjes_bandwidth_eta: 0.1,
            },
        );
        let eg = ql_portfolio::evaluate_tsay_egarch_and_hamilton_markov_switching(
            &ql_portfolio::TsayEgarchHamiltonRegimeSpec {
                returns: vec![0.01 + r * 0.01, -0.02 + credit_bump, 0.015 + vols * 0.001],
                initial_variance: (0.04 + vols * 0.001).max(0.001),
                egarch_omega: -0.1,
                egarch_alpha: 0.1,
                egarch_leverage_gamma: -0.1,
                egarch_beta: 0.9,
                calm_mean: 0.0,
                calm_vol: 0.1,
                crisis_mean: 0.0,
                crisis_vol: 0.3,
                stay_calm_p00: 0.95,
                stay_crisis_p11: 0.9,
            },
        );
        let iv = ql_dml::evaluate_wooldridge_greene_2sls_iv_and_did(
            &ql_dml::WooldridgeGreene2slsDidSpec {
                endogenous_x: vec![1.0 + credit_bump, 2.0 + r * 0.01, 3.0 + vols * 0.01],
                instrument_z: vec![1.0, 2.1, 2.9],
                outcome_y: vec![1.2, 2.1, 3.2],
                treated_pre_mean: 1.0,
                treated_post_mean: 1.2,
                control_pre_mean: 1.0,
                control_post_mean: 1.1,
            },
        );
        let mc = ql_portfolio::evaluate_mcneil_frey_copula_tail_and_euler_allocation(
            &ql_portfolio::McneilFreyCopulaEulerAllocationSpec {
                clayton_theta: 2.0,
                gumbel_theta: 1.5,
                student_t_df_nu: 5.0,
                student_t_corr_rho: 0.3,
                desk_weights: [0.25 + credit_bump * 0.1, 0.25, 0.25, 0.25],
                desk_expected_losses: [0.01; 4],
                loss_covariance_4x4: [
                    [0.04, 0.01, 0.0, 0.0],
                    [0.01, 0.05, 0.01, 0.0],
                    [0.0, 0.01, 0.06, 0.01],
                    [0.0, 0.0, 0.01, 0.03],
                ],
                es_tail_multiplier: 2.0,
            },
        );
        let om = ql_portfolio::evaluate_bacon_omega_sortino_and_burke_drawdown(
            &ql_portfolio::BaconOmegaSortinoBurkeDrawdownSpec {
                period_returns: vec![0.01 + r * 0.01, -0.005 + credit_bump, 0.02 + vols * 0.001],
                mar_threshold: 0.0,
                risk_free_rate: 0.002,
            },
        );
        let rd = ql_models::evaluate_damodaran_rd_capitalization_and_patent_real_option(
            &ql_models::DamodaranRdCapitalizationRealOptionSpec {
                reported_ebit: 10.0,
                reported_invested_capital: 100.0,
                effective_tax_rate: 0.25,
                current_rd_expense: 5.0,
                historical_rd_expenses: vec![4.0, 3.0, 2.0],
                patent_project_pv: 20.0,
                patent_development_cost: 10.0,
                patent_life_years: 5.0,
                project_volatility: 0.3,
                risk_free_rate: r.max(0.001),
                abandonment_salvage_value: 2.0,
            },
        );
        let ad = ql_models::evaluate_mascolell_cochrane_arrow_debreu_sdf(
            &ql_models::MasColellCochraneArrowDebreuSdfSpec {
                state_probabilities: [0.3, 0.4, 0.3],
                consumption_growth_ratios: [0.98, 1.02 + credit_bump, 1.05 + r * 0.01],
                equity_payoff_states: [0.8, 1.1, 1.4],
                subjective_discount_beta: 0.99,
                crra_risk_aversion_gamma: 2.0,
            },
        );
        c.forward_swap_rate
            + cr.structural_credit_spread_bps * 1e-4
            + ca.expected_gamma_scalp_pnl
            + sp.specified_pool_payup_points * 1e-2
            + og.inventory_skewed_reservation_mid * 1e-4
            + rie.rie_condition_number * 1e-2
            + eg.terminal_filtered_crisis_probability
            + iv.two_stage_ls_iv_beta * 1e-2
            + mc.total_portfolio_es
            + om.keating_shadwick_omega_ratio * 1e-2
            + rd.total_strategic_expanded_npv * 1e-2
            + ad.equilibrium_equity_price * 1e-2
    }
    pub fn rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric(&self) -> f64 {
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else if self.scenario == ScenarioKind::StagflationPlus150Bp {
            0.005
        } else {
            0.0
        };
        Self::rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_live_metric(&self) -> f64 {
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else if self.scenario == ScenarioKind::StagflationPlus150Bp {
            0.005
        } else {
            0.0
        };
        Self::cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump;
        let inf = m.inflation_quotes.first().map(|(_, v)| *v).unwrap_or(0.02);
        let k = ql_curves::evaluate_fabozzi_tuckman_key_rate_and_lsc_immunization(
            &ql_curves::FabozziKeyRateLscImmunizationSpec {
                cashflow_tenors: vec![1.0, 3.0, 5.0, 10.0],
                cashflows: vec![2.0, 2.0, 2.0, 102.0],
                zero_rates: vec![r, r + 0.002, r + 0.004, r + 0.006],
                key_rate_tenors: vec![1.0, 5.0, 10.0],
                ns_decay_lambda: 0.5 + inf,
                liability_horizon_years: 5.0,
                barbell_short_tenor: 1.0,
                barbell_long_tenor: 10.0,
                twist_slope_shock_bps: 25.0,
                butterfly_curvature_shock_bps: 15.0,
            },
        );
        let cmo = ql_fixed_income::evaluate_davidson_levin_cmo_pac_companion_waterfall(
            &ql_fixed_income::DavidsonLevinCmoPacCompanionSpec {
                initial_pool_balance: 100.0,
                gross_wac_rate: r + 0.03,
                net_coupon_rate: r + 0.02,
                num_periods: 24,
                pac_lower_psa: 100.0,
                pac_upper_psa: 300.0,
                realized_psa: 180.0 + credit_bump * 100.0,
                discount_yield: r,
            },
        );
        let cb = ql_fixed_income::evaluate_asquith_sirbu_shreve_cb_call_delay_and_dynkin_game(
            &ql_fixed_income::AsquithSirbuShreveCbGameSpec {
                firm_value: 130.0,
                stock_price: 100.0,
                conversion_ratio: 1.0,
                annual_coupon_amount: 3.0,
                annual_dividend_per_share: 1.0,
                corporate_tax_rate: 0.25,
                effective_call_price: 105.0,
                call_notice_period_years: 30.0 / 365.0,
                firm_volatility: m.eq_vol.max(0.05),
                risk_free_rate: r,
                firm_payout_rate: 0.02,
                dilution_gamma: 0.8,
            },
        );
        let asu = ql_fixed_income::evaluate_fasb_asu2024_induced_conversion_and_cpdi_oid(
            &ql_fixed_income::FasbAsu2024InducedConversionCpdiSpec {
                net_income: 12.0,
                basic_shares: 100.0,
                bond_par: 100.0,
                stated_coupon_rate: 0.03,
                comparable_noncontingent_yield: r + 0.01,
                tax_rate: 0.25,
                original_conversion_shares: 1.0,
                induced_sweetener_extra_shares: 0.05,
                induced_sweetener_cash: 1.0,
                stock_price_at_acceptance: m.eq_spot,
                maturity_years: 5.0,
            },
        );
        let modif = ql_credit::evaluate_rsm_asc470_debt_modification_and_restructuring(
            &ql_credit::RsmAsc470DebtModificationSpec {
                carrying_amount: 100.0,
                original_par: 100.0,
                original_coupon_rate: 0.03,
                original_effective_rate: r,
                modified_par: 100.0,
                modified_coupon_rate: 0.025,
                remaining_years: 5.0,
                lender_fee_paid: 1.0,
                conversion_option_fv_change: 0.5,
                market_exit_yield: r + 0.01,
            },
        );
        let colva = ql_xva::evaluate_gregory_csa_threshold_mta_and_ctd_colva(
            &ql_xva::GregoryCsaThresholdMtaCtdColvaSpec {
                time_grid: vec![1.0, 2.0, 3.0, 5.0],
                expected_mtm_profile: vec![5.0, 8.0, 10.0, 12.0],
                exposure_vol: 3.0 + m.eq_vol,
                mpor_years: 10.0 / 252.0,
                csa_threshold: 2.0,
                minimum_transfer_amount: 1.0,
                independent_amount: 0.5,
                hazard_rate: 0.02 + credit_bump,
                lgd: 0.6,
                risk_free_rate: r,
                ccy1_collateral_spread: 0.01,
                ccy2_collateral_spread: 0.008,
                collateral_spread_vol: 0.01,
            },
        );
        let evt = ql_portfolio::evaluate_embrechts_hill_pickands_and_subexponential_tail(
            &ql_portfolio::EmbrechtsHillPickandsSubexponentialSpec {
                loss_samples: vec![
                    160.0, 110.0, 82.0, 66.0, 56.0, 50.0, 46.0, 43.0, 40.0, 38.0, 36.0, 35.0,
                ],
                tail_order_k: 2,
                extreme_tail_prob: 0.01,
                portfolio_positions_count: 8,
            },
        );
        let tlf = ql_portfolio::evaluate_bouchaud_potters_power_law_cumulant_and_smile_scaling(
            &ql_portfolio::BouchaudPottersCumulantSmileSpec {
                single_period_vol: m.eq_vol,
                single_period_skewness: -0.2,
                single_period_excess_kurtosis: 3.0,
                vol_of_vol_kurtosis_component: 2.0,
                vol_decorrelation_steps: 10.0,
                horizon_steps_n: 20.0,
                log_moneyness: 0.0,
                confidence_z: 2.33,
            },
        );
        let attr = ql_portfolio::evaluate_karnosky_singer_multicurrency_and_menchero_attribution(
            &ql_portfolio::KarnoskySingerMencheroAttributionSpec {
                portfolio_asset_weights: vec![0.6, 0.4],
                benchmark_asset_weights: vec![0.5, 0.5],
                portfolio_currency_weights: vec![0.6, 0.4],
                benchmark_currency_weights: vec![0.5, 0.5],
                local_portfolio_returns: vec![0.08, 0.03],
                local_benchmark_returns: vec![0.06, 0.04],
                local_eurodollar_rates: vec![0.02, 0.025],
                base_eurodollar_rate: r,
                spot_fx_returns: vec![0.01, -0.005],
                multi_period_port_returns: vec![0.02, 0.03],
                multi_period_bench_returns: vec![0.018, 0.025],
            },
        );
        let tm = ql_portfolio::evaluate_dgtw_characteristic_and_treynor_mazuy_timing(
            &ql_portfolio::DgtwTreynorMazuyTimingSpec {
                current_weights: vec![0.6, 0.4],
                lagged_weights: vec![0.5, 0.5],
                stock_returns: vec![0.08, 0.03],
                current_char_bench_returns: vec![0.06, 0.04],
                lagged_char_bench_returns: vec![0.05, 0.04],
                portfolio_excess_series: vec![0.01, 0.02, 0.04, 0.08],
                market_excess_series: vec![0.0, 0.02, 0.04, 0.08],
            },
        );
        let bh = ql_portfolio::evaluate_brock_hommes_adaptive_belief_switching(
            &ql_portfolio::BrockHommesAdaptiveBeliefSpec {
                initial_price_deviation: 0.1,
                risk_free_rate: r,
                fundamentalist_damping_v1: 0.8,
                chartist_trend_extrapolation_g: 1.2,
                fundamentalist_info_cost_c1: 0.1,
                intensity_of_choice_beta: 8.0,
                fitness_memory_eta: 0.8,
                risk_aversion_sigma2: 1.0,
                shock_innovations: vec![0.01, -0.005, 0.003],
            },
        );
        let pen = ql_models::evaluate_penman_reoi_and_ohlson_juettner_aeg_valuation(
            &ql_models::PenmanReoiOhlsonJuettnerSpec {
                net_operating_assets_0: 100.0,
                net_financial_obligations_0: 30.0,
                operating_income_after_tax_1: 12.0,
                operating_income_after_tax_2: 13.0,
                net_financial_expense_after_tax_1: 2.0,
                wacc: r + 0.02,
                terminal_reoi_growth: 0.02,
                current_share_price: m.eq_spot,
                forward_eps_1: 5.0,
                forward_eps_2: 5.5,
                forward_dps_1: 2.0,
                perpetual_aeg_growth: 0.03,
            },
        );
        k.bond_present_value
            + cmo.pac_present_value
            + cb.perpetual_dynkin_cb_fair_value
            + cb.asquith_safety_call_trigger_stock_price * 1e-2
            + cb.conversion_value_to_call_ratio
            + asu.basic_eps
            + asu.asu2020_06_instrument_c_net_share_eps
            + asu.asu2024_04_induced_conversion_expense * 1e-2
            + asu.cpdi_projected_contingent_terminal_payment * 1e-3
            + modif.pv_modified_cashflows_at_orig_eir
            + modif.prospective_modified_eir
            + colva.ctd_multi_currency_colva
            + evt.weissman_extreme_var * 1e-3
            + tlf.bouchaud_potters_non_gaussian_var
            + attr.cumulative_active_return
            + attr.karnosky_singer_security_selection
            + attr.karnosky_singer_currency_management
            + attr.menchero_natural_scaling_multiplier * 1e-3
            + tm.treynor_mazuy_timing_gamma
            + tm.dgtw_characteristic_selectivity_cs
            + bh.terminal_price_deviation
            + pen.reoi_intrinsic_equity_value * 1e-2
            + pen.ohlson_juettner_implied_cost_of_equity
    }

    pub fn fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_live_metric(&self) -> f64 {
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_live_metric_from_market(
            &self.scenario_market(),
            bump,
        )
    }

    pub fn veronesi_conv_lando_gen_gatheral_corridor_live_metric_from_market(
        m: &MarketSnapshot,
    ) -> f64 {
        Self::veronesi_conv_lando_gen_gatheral_corridor_live_metric_from_market_with_bump(m, 0.0)
    }

    pub fn veronesi_conv_lando_gen_gatheral_corridor_live_metric_from_market_with_bump(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03);
        let spot = m.eq_spot.max(1e-6);
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let a = ql_fixed_income::evaluate_veronesi_fabozzi_futures_convexity_and_prepayment(
            &ql_fixed_income::VeronesiFabozziConvexityPrepaymentSpec {
                forward_rate: r,
                short_rate_vol: vol,
                mean_reversion_kappa: 0.1,
                fixing_time_t1: 1.0,
                maturity_time_t2: 2.0,
                wac_coupon: 0.05,
                current_mortgage_rate: r + inf,
                min_cpr_turnover: 0.02,
                max_cpr_refi: 0.3,
                refi_steepness_gamma: 8.0,
                refi_incentive_threshold: 0.0,
                pool_age_months: 24.0,
                psa_ramp_months: 30.0,
            },
        );
        let b = ql_credit::compute_lando_israel_generator_logarithm_and_migration(
            &ql_credit::LandoIsraelGeneratorLogarithmSpec {
                transition_matrix_3x3: [[0.9, 0.08, 0.02], [0.03, 0.92, 0.05], [0.0, 0.0, 1.0]],
                series_order_k: 12,
                fractional_horizon_years: 0.5,
                recovery_rate: 0.4,
                risk_free_rate: r + credit_bump,
            },
        );
        let c = ql_vol::evaluate_gatheral_demeterfi_corridor_variance_and_skew(
            &ql_vol::GatheralDemeterfiCorridorVarianceSkewSpec {
                forward: spot,
                maturity: 1.0,
                atm_vol: vol,
                atm_skew: -0.2,
                smile_curvature: 0.1,
                corridor_lower_strike: spot * 0.9,
                corridor_upper_strike: spot * 1.1,
                sampling_steps_n: 12,
                vol_of_vol_eta: vol,
                lsv_mixing_weight: 0.3,
            },
        );
        let d = ql_exotics::price_shepp_shiryaev_shreve_perpetual_russian_and_mckean(
            &ql_exotics::SheppShiryaevShrevePerpetualRussianSpec {
                spot,
                running_max: spot * 1.1,
                fixed_put_strike: spot,
                risk_free_rate: r,
                dividend_yield: inf,
                volatility: vol,
            },
        );
        let e = ql_levy::evaluate_cont_tankov_bates_svj_characteristic_and_pide(
            &ql_levy::ContTankovBatesSvjPideSpec {
                spot,
                strike: spot,
                maturity: 1.0,
                risk_free_rate: r,
                dividend_yield: inf,
                v0: vol * vol,
                kappa: 1.0,
                theta: vol * vol,
                sigma_v: vol,
                rho: -0.5,
                jump_intensity_lambda: 0.2,
                jump_mean_mu_j: -0.1,
                jump_std_sigma_j: 0.2,
            },
        );
        let f = ql_fixed_income::evaluate_koziol_asquith_fasb_strategic_convertible(
            &ql_fixed_income::KoziolAsquithFasbStrategicConvertibleSpec {
                firm_value: spot * 100.0,
                debt_face_value: spot * 80.0,
                existing_shares: 100.0,
                new_shares_blockholder1: 10.0,
                new_shares_blockholder2: 8.0,
                coupon_rate: 0.04,
                corporate_tax_rate: 0.25,
                dividend_yield: inf,
                stock_price: spot,
                base_conversion_ratio: 1.0,
                induced_sweetened_conversion_ratio: 1.1,
                cash_sweetener_amount: 1.0,
                call_price: spot * 1.2,
                firm_volatility: vol,
                risk_free_rate: r,
                maturity: 5.0,
            },
        );
        let g = ql_models::solve_carmona_delarue_mean_field_systemic_interbank_game(
            &ql_models::CarmonaDelarueMeanFieldSystemicSpec {
                initial_log_reserve_spread: 0.01 + credit_bump,
                interbank_mean_reversion_a: 0.4,
                control_quadratic_penalty_q: 0.2,
                systemic_risk_preference_epsilon: 0.3,
                terminal_systemic_penalty_c: 1.0,
                idiosyncratic_vol_sigma: vol,
                common_noise_corr_rho0: 0.2,
                distress_barrier_d: -0.1,
                horizon_t: 1.0,
            },
        );
        let h = ql_xva::evaluate_bcbs_d424_output_floor_and_sma_op_risk(
            &ql_xva::BcbsD424OutputFloorAndSmaOpRiskSpec {
                irb_credit_rwa: 100.0,
                ima_market_rwa: 50.0,
                standardized_credit_rwa: 120.0,
                standardized_market_rwa: 60.0,
                standardized_cva_rwa: 20.0,
                output_floor_multiplier: 0.725,
                interest_lease_dividend_component: 20.0,
                services_component: 15.0,
                financial_component: 30.0,
                average_annual_operational_loss: 2.0,
                cet1_capital: 20.0,
                at1_and_tier2_capital: 10.0,
                eligible_tlac_senior_debt: 40.0,
            },
        );
        let i = ql_portfolio::solve_cartea_jaimungal_penalva_lob_fill_hjb(
            &ql_portfolio::CarteaJaimungalPenalvaLobFillHjbSpec {
                initial_inventory: 10,
                horizon_seconds: 60.0,
                base_arrival_rate_lambda: 0.5,
                fill_decay_kappa: 1.0,
                running_inventory_penalty_phi: 0.01,
                terminal_liquidation_penalty_alpha: 0.1,
                ofi_short_term_drift: 0.02,
                adverse_selection_cost: 0.01,
            },
        );
        let j = ql_portfolio::evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha(
            &ql_portfolio::MadhavanAbnerHillEtfApArbitrageSpec {
                etf_market_price: spot * 1.01,
                reported_inav: spot,
                futures_beta_return: 0.01,
                stale_underlying_weight: 0.4,
                ap_creation_fee_bps: 5.0,
                basket_bid_ask_spread_bps: 8.0,
                creation_unit_volatility_bps: 10.0,
                ap_risk_aversion: 0.5,
                embedded_unrealized_gain_fraction: 0.3,
                annual_redemption_turnover: 0.2,
                capital_gains_tax_rate: 0.25,
            },
        );
        let k = ql_portfolio::evaluate_grinold_qian_christopherson_ic_decay_and_timing(
            &ql_portfolio::GrinoldQianChristophersonAlphaDecaySpec {
                fast_signal_ic0: 0.1,
                fast_signal_half_life_months: 3.0,
                slow_signal_ic0: 0.08,
                slow_signal_half_life_months: 12.0,
                signal_correlation: 0.2,
                rebalance_horizon_months: 1.0,
                transaction_cost_penalty: 0.01,
                factor_active_exposures: vec![0.1, -0.05],
                factor_benchmark_returns: vec![0.02, 0.01],
                normal_factor_exposures: vec![0.05, -0.02],
            },
        );
        let l = ql_dml::evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad(
            &ql_dml::GlassermanSaddlepointIsStratifiedLhsSpec {
                spot,
                deep_otm_strike: spot * 1.3,
                risk_free_rate: r,
                volatility: vol,
                maturity: 1.0,
                num_strata: 16,
            },
        );
        a.implied_futures_rate
            + b.ig_to_default_fractional_pd
            + c.full_variance_swap_fair_vol
            + d.perpetual_russian_option_value
            + e.bates_svj_call_price
            + f.strategic_convertible_fair_value * 1e-3
            + g.mean_field_value_function_cost
            + h.floored_cet1_ratio
            + i.expected_terminal_inventory * 1e-3
            + j.net_ap_creation_redemption_profit_bps * 1e-3
            + k.total_factor_contribution
            + l.importance_sampled_call_price
            + inf
    }

    pub fn veronesi_conv_lando_gen_gatheral_corridor_live_metric(&self) -> f64 {
        let credit_bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::veronesi_conv_lando_gen_gatheral_corridor_live_metric_from_market_with_bump(
            &self.scenario_market(),
            credit_bump,
        )
    }

    pub fn bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_live_metric_from_market(
        m: &MarketSnapshot,
    ) -> f64 {
        Self::bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_live_metric_from_market_with_bump(m, 0.0)
    }

    fn bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_live_metric_from_market_with_bump(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let r = m.ois_quotes.first().copied().unwrap_or(0.03) + credit_bump;
        let spot = m.eq_spot;
        let vol = m.eq_vol.max(1e-4);
        let credit = (m.fwd_quotes.first().copied().unwrap_or(r) - r)
            .abs()
            .max(0.01);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let irr = evaluate_bcbs_d368_six_scenario_irrbb(&BcbsD368SixScenarioIrrbbSpec {
            bucket_midpoints_years: vec![0.5, 2.0, 5.0],
            asset_repricing_cashflows: vec![10.0 + inf, 8.0, 5.0],
            liability_repricing_cashflows: vec![7.0, 6.0, 4.0],
            baseline_zero_rates: vec![r, r + 0.002, r + 0.004],
            parallel_shock_bps: 200.0,
            short_shock_bps: 250.0,
            long_shock_bps: 100.0,
            decay_tau_years: 4.0,
            tier1_capital: 100.0,
            pass_through_deposit_beta: 0.5,
        });
        let rfr = evaluate_turfus_willems_rfr_in_arrears_caplet_smile(
            &TurfusWillemsRfrInArrearsSmileSpec {
                accrual_start_t1: 0.1,
                accrual_end_t2: 0.6,
                forward_rfr_rate: r,
                strike_rate: r * 0.9,
                discount_factor_t2: (-r * 0.6).exp(),
                sabr_alpha: vol * 0.1,
                sabr_beta: 0.5,
                sabr_rho: -0.2,
                sabr_nu: vol,
                hull_white_mean_reversion: 0.1,
                notional: 1e6,
            },
        );
        let fx_spot = 1.1 + spot * 1e-4;
        let xccy =
            evaluate_baran_bianchetti_mtm_xccy_basis_swap(&BaranBianchettiMtmXccyBasisSpec {
                spot_fx_dom_per_for: fx_spot,
                foreign_notional: 100.0,
                payment_times: vec![0.5, 1.0, 1.5],
                dom_ois_discount_factors: vec![(-r * 0.5).exp(), (-r).exp(), (-r * 1.5).exp()],
                for_xccy_discount_factors: vec![0.98, 0.96, 0.94],
                dom_forward_rates: vec![r; 3],
                for_forward_rates: vec![r - 0.005; 3],
                quoted_dom_basis_spread: credit * 0.01,
                dom_rate_vol: vol,
                fx_spot_vol: vol,
                corr_dom_rate_fx: -0.3,
            });
        let bg = evaluate_bourgey_gatheral_bergomi_guyon_expansion(
            &BourgeyGatheralBergomiGuyonSmileSpec {
                forward_price: spot,
                strike_price: spot,
                maturity: 1.0,
                initial_instantaneous_vol: vol,
                vol_of_vol_eta: vol,
                spot_vol_corr_rho: -0.5,
                hurst_exponent: 0.1,
                mean_reversion_lambda: 1.0,
            },
        );
        let coco = evaluate_feds_despiegeleer_at1_coco_dual_trigger(
            &FedsDeSpiegeleerAt1CocoDualTriggerSpec {
                notional: 100.0,
                coupon_rate: 0.08,
                maturity: 5.0,
                risk_free_rate: r,
                initial_cet1_ratio: 0.12 + credit,
                mechanical_cet1_trigger: 0.07,
                ponv_cet1_trigger: 0.05,
                supervisory_ponv_hazard_lambda: credit,
                asset_volatility: vol,
                cet1_asset_elasticity: 1.0,
                initial_stock_price: spot,
                conversion_price: spot * 0.8,
                floor_stock_price: spot * 0.2,
                existing_common_shares: 100.0,
                write_down_fraction: 0.5,
                is_equity_conversion: true,
            },
        );
        let xv = evaluate_albanese_crepey_zeron_balance_sheet_xva(
            &AlbaneseCrepeyZeronBalanceSheetXvaSpec {
                maturity: 2.0,
                num_time_steps: 8,
                risk_free_rate: r,
                cpty_hazard_rate: credit,
                cpty_recovery: 0.4,
                bank_hazard_rate: credit * 0.5,
                bank_recovery: 0.4,
                bank_liquidity_basis: credit * 0.2,
                im_funding_spread: credit,
                shareholder_hurdle_rate: 0.1,
                base_unmargined_epe: 10.0 + inf,
                base_ene: 4.0,
                simm_delta_weight: 0.6,
                simm_vega_weight: 0.4,
                simm_delta_vega_corr: 0.2,
                regulatory_capital_ratio: 0.08,
            },
        );
        let rb =
            evaluate_jacquier_chebyshev_rough_vol_surrogate(&JacquierDeepRoughVolSurrogateSpec {
                forward: spot,
                maturity: 1.0,
                initial_var_xi0: vol * vol,
                hurst_h: 0.1,
                vol_of_vol_eta: vol,
                spot_vol_corr_rho: -0.7,
                target_market_vols_3_strikes: [vol * 1.1, vol, vol * 1.05],
                log_moneyness_strikes: [-0.1, 0.0, 0.1],
            });
        let sp = evaluate_choi_li_deng_multi_asset_spread_basket(
            &ChoiLiDengMultiAssetSpreadBasketSpec {
                long_leg_forward: spot,
                long_leg_weight: 1.0,
                long_leg_vol: vol,
                short_leg_forwards: vec![spot * 0.95],
                short_leg_weights: vec![1.0],
                short_leg_vols: vec![vol * 1.1],
                corr_long_short: vec![0.4],
                corr_short_short_equicorr: 0.0,
                strike: spot * 0.02,
                maturity: 1.0,
                risk_free_rate: r,
            },
        );
        let bs = evaluate_bess_analytical_sdp_degradation_storage(
            &BessAnalyticalSdpDegradationStorageSpec {
                energy_capacity_mwh: 100.0,
                max_power_mw: 50.0,
                initial_soc_fraction: 0.5,
                min_soc_fraction: 0.1,
                max_soc_fraction: 0.9,
                charge_efficiency: 0.9,
                discharge_efficiency: 0.9,
                battery_capex_usd: 100000.0,
                rated_cycle_life: 4000.0,
                peak_mean_price_mwh: spot,
                offpeak_mean_price_mwh: spot * 0.5,
                price_volatility_mwh: spot * 0.2,
                horizon_days: 365,
                discount_rate_annual: r,
            },
        );
        let dk = evaluate_mitchell_pedersen_slow_capital_and_kratz_dark_pool(
            &MitchellPedersenKratzDarkPoolSpiralSpec {
                initial_inventory_shares: 1e5,
                reference_price: spot,
                execution_horizon_days: 5.0,
                daily_volatility: vol,
                lit_temporary_impact_eta: 1e-6,
                lit_permanent_impact_gamma: 2e-6,
                risk_aversion_lambda: 0.1,
                dark_pool_poisson_intensity: 0.3,
                dark_pool_adverse_selection_bps: 5.0,
                specialist_immediate_capital: 1e7,
                slow_moving_external_capital: 2e7,
                capital_arrival_speed_beta: 0.4,
            },
        );
        let rm = evaluate_follmer_schied_convex_risk_and_diagnostic_expectations(
            &FollmerSchiedGennaioliDiagnosticRiskSpec {
                scenario_pnl: vec![spot * 0.01, -spot * 0.02, -spot * 0.05],
                baseline_probabilities: vec![0.5, 0.3, 0.2],
                entropic_ambiguity_theta: 2.0,
                kl_divergence_budget_eta: 0.1,
                diagnostic_overreaction_theta: 1.5,
                boom_news_drift_shock: inf,
                bust_news_drift_shock: credit,
            },
        );
        let pc =
            evaluate_huge_savine_differential_pca_and_twin_dml(&HugeSavineDifferentialPcaSpec {
                state_samples: vec![vec![0.2, 2.0], vec![0.3, -1.8], vec![-0.25, 1.5]],
                payoff_samples: vec![spot, spot * 0.98, spot * 1.02],
                aad_gradient_samples: vec![vec![2.2, 0.1], vec![1.8, -0.15], vec![-2.0, 0.12]],
                twin_fd_gradient_samples: vec![vec![2.0, 0.12], vec![1.7, -0.1], vec![-1.9, 0.1]],
                retained_rank: 1,
                ridge_regularization: 1e-6,
            });
        irr.baseline_eve
            + rfr.in_arrears_caplet_pv * 1e-6
            + xccy.mtm_swap_npv_in_dom
            + bg.atm_implied_vol_second_order
            + coco.at1_coco_fair_value * 1e-3
            + xv.ftd_cva
            + rb.surrogate_smile_vols[0]
            + sp.choi_gauss_hermite_price * 1e-2
            + bs.total_bess_npv_usd * 1e-5
            + dk.kratz_schoneborn_dark_lit_cost * 1e-3
            + rm.follmer_schied_entropic_risk_measure
            + pc.truncation_gradient_mse_bound
            + inf
    }

    pub fn bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_live_metric(&self) -> f64 {
        let credit_bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_live_metric_from_market_with_bump(
            &self.scenario_market(),
            credit_bump,
        )
    }

    pub fn cochrane_ftpl_g2_twist_essvi_offshore_cb_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let credit_bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::cochrane_ftpl_g2_twist_essvi_offshore_cb_live_metric_from_market(&market, credit_bump)
    }

    pub fn mf_rates_wildcard_ctd_pivot_tarf_bergomi_live_metric_from_market(
        m: &MarketSnapshot,
        credit_bump: f64,
    ) -> f64 {
        let r = m.ois_quotes[0] + credit_bump;
        let spot = m.eq_spot;
        let vol = m.eq_vol.max(1e-4);
        let inf = m.inflation_quotes.first().map(|(_, z)| *z).unwrap_or(0.023);
        let times = vec![0.25, 0.5, 0.75, 1.0];
        let mf = evaluate_hunt_kennedy_markov_functional_range_accrual(
            &HuntKennedyMarkovFunctionalSpec {
                fixing_times: times.clone(),
                initial_forwards: vec![r; 4],
                market_caplet_vols: vec![vol; 4],
                smile_curvature_eta: 0.1 * inf,
                mean_reversion: 0.03,
                range_lower: -0.01,
                range_upper: 0.03,
                accrual_coupon_rate: 0.02 + inf,
                call_strike_spread: 0.005,
            },
        );
        let wc = evaluate_tuckman_burghardt_treasury_wildcard_and_ctd_switch(
            &TuckmanBurghardtTreasuryWildcardBasisSpec {
                futures_settlement_price: 100.0,
                delivery_horizon_years: 0.25,
                general_collateral_repo_rate: r,
                repo_specialness_spread: 0.002,
                bond1_clean_price: 99.0,
                bond1_conversion_factor: 0.85,
                bond1_modified_duration: 6.0,
                bond1_yield_beta: 1.0,
                bond2_clean_price: 101.0,
                bond2_conversion_factor: 0.9,
                bond2_modified_duration: 7.0,
                bond2_yield_beta: 1.1,
                afternoon_wildcard_vol_annual: vol,
                eom_delivery_window_days: 5,
            },
        );
        let tarf = evaluate_wystup_pivot_tarf_and_corridor_strip(&WystupPivotTarfCorridorSpec {
            spot,
            strike_low: 0.95 * spot,
            strike_pivot: spot,
            strike_high: 1.05 * spot,
            downside_gearing: 1.2,
            upside_gearing: 1.0,
            target_redemption_cap: 0.05,
            domestic_rate: r,
            foreign_rate: 0.01,
            volatility: vol,
            fixing_times: times.clone(),
        });
        let berg = evaluate_bergomi_forward_start_smile_and_napoleon_convexity(
            &BergomiForwardStartNapoleonVolSpec {
                initial_forward_variance: vol * vol,
                vol_of_vol_omega: 0.3,
                kappa_short_k1: 3.0,
                kappa_long_k2: 0.5,
                long_factor_weight_theta: 0.5,
                factor_corr_rho12: 0.3,
                spot_vol_corr_rho_s1: -0.5,
                spot_vol_corr_rho_s2: -0.3,
                forward_start_t1: 0.5,
                forward_end_t2: 1.0,
                forward_moneyness_k: 1.0,
                napoleon_base_coupon: 0.01,
                napoleon_global_floor: 0.0,
                napoleon_subperiods_m: 4,
            },
        );
        let sch = evaluate_schonbucher_defaultable_range_accrual_and_spread_cap(
            &SchonbucherDefaultableRangeAccrualCapSpec {
                notional: 100.0,
                payment_times: times.clone(),
                risk_free_rate: r,
                initial_hazard_rates: vec![0.02 + credit_bump; 4],
                base_recovery_rate: 0.4,
                spread_volatilities: vec![vol; 4],
                spread_cap_strike: 0.03,
                accrual_spread_lower: 0.0,
                accrual_spread_upper: 0.05,
                range_accrual_coupon_rate: 0.02,
                recovery_spread_correlation: -0.2,
            },
        );
        let cva = evaluate_green_frtb_cva_hedge_relief_and_kva(&GreenFrtbCvaHedgeCapitalSpec {
            counterparty_sacr_eads: vec![100.0 * (1.0 + credit_bump), 80.0 * (1.0 + vol)],
            counterparty_effective_maturities: vec![3.0 + inf, 5.0 + vol],
            counterparty_risk_weights: vec![0.02 + credit_bump, 0.03 + inf.abs()],
            single_name_hedge_notionals: vec![20.0 * (1.0 + spot.abs() * 1e-3)],
            single_name_hedge_maturities: vec![3.0],
            hedge_reference_correlations_rhc: vec![0.8],
            index_hedge_notional: 10.0 * (1.0 + r.abs()),
            index_hedge_maturity: 5.0,
            index_risk_weight: 0.02,
            sa_cva_spread_deltas: vec![1.0, 2.0],
            sa_cva_vega_sensitivities: vec![0.5],
            hurdle_rate_gamma_k: 0.1,
            hedge_margin_funding_spread: 0.01 + vol,
            discount_rate: r.max(1e-4),
        });
        let lb = evaluate_heynen_kat_partial_fractional_lookback_and_extreme_spread(
            &HeynenKatPartialFractionalLookbackSpec {
                spot,
                running_min: 0.9 * spot,
                running_max: 1.1 * spot,
                fractional_multiplier_lambda: 1.2,
                monitoring_end_t1: 0.5,
                maturity_t2: 1.0,
                fixed_strike: spot,
                rate: r,
                dividend_yield: 0.0,
                volatility: vol,
            },
        );
        let vta = evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling(
            &EydelandWolyniecRampOutageVtaSpec {
                power_forwards: vec![60.0, 62.0],
                gas_forwards: vec![30.0, 31.0],
                spinning_reserve_prices: vec![3.0, 3.0],
                dispatch_times: times[..2].to_vec(),
                capacity_mw: 100.0,
                block_hours: 1.0,
                full_load_heat_rate: 2.0,
                part_load_heat_rate: 2.2,
                min_stable_load_fraction: 0.5,
                vom_cost_per_mwh: 1.0,
                startup_cost: 2.0,
                fast_ramp_penalty_cost: 1.0,
                outage_failure_rate_lambda: 0.05,
                outage_repair_rate_mu: 0.95,
                spread_volatility: vol,
                risk_free_rate: r,
            },
        );
        let fl =
            evaluate_grinold_kahn_clarke_fundamental_law(&GrinoldKahnClarkeFundamentalLawSpec {
                standardized_scores_z: vec![1.0 + inf, -0.5 + credit_bump, 0.7 + spot * 1e-3],
                residual_volatilities: vec![
                    vol,
                    (vol + r.abs()).max(1e-4),
                    (vol + inf.abs()).max(1e-4),
                ],
                benchmark_weights: vec![0.3, 0.4, 0.3],
                information_coefficient_ic: 0.1,
                rebalances_per_year: 12,
                signal_autocorrelation_rho: 0.0,
                target_tracking_error: 0.05,
                max_active_weight_bound: 0.1,
            });
        let acd = evaluate_engle_russell_acd_trade_duration_and_hazard(
            &EngleRussellAcdTradeDurationSpec {
                inter_trade_durations: vec![1.0 + vol, 1.2 + r.abs(), 0.9 + spot.abs() * 1e-3],
                omega: 0.2 + r.abs(),
                alpha: 0.2 + inf.abs(),
                beta: 0.5 + credit_bump.abs(),
                weibull_shape_gamma: 1.0 + vol,
                initial_expected_duration: 1.0,
                elapsed_time_since_last_trade: 0.5 + inf.abs(),
                base_price_impact_bps: 2.0 + spot.abs() * 1e-3,
                duration_impact_elasticity: 0.5,
            },
        );
        let hlz = evaluate_harvey_liu_zhu_multiple_testing_and_mintrl(
            &HarveyLiuZhuMultipleTestingSharpeSpec {
                focal_sharpe_ratio: 0.8 + inf + credit_bump,
                benchmark_sharpe_sr0: 0.0,
                sample_length_t: 252,
                return_skewness: -0.5,
                return_kurtosis_pearson: 4.0,
                num_independent_trials_m: 20,
                cross_sectional_sharpe_std: 0.3,
                confidence_z_alpha: 1.96 + vol,
                candidate_t_stats: vec![1.5 + spot * 1e-3, 2.0 + vol],
            },
        );
        let ml = evaluate_giles_capriotti_multilevel_mc_aad(&GilesCapriottiMultilevelMcAadSpec {
            spot,
            strike: spot,
            rate: r,
            volatility: vol,
            local_vol_skew: -0.2,
            maturity: 1.0,
            base_steps_m0: 8,
            target_rmse_epsilon: 0.01,
        });
        mf.non_callable_range_accrual_npv.val()
            + wc.theoretical_net_basis_ctd.val() * 1e-2
            + tarf.pivot_tarf_client_npv.val()
            + berg.bergomi_stoch_vol_napoleon_coupon.val()
            + sch.defaultable_range_accrual_note_npv.val() * 1e-2
            + cva.ba_cva_total_capital_charge.val() * 1e-2
            + lb.heynen_kat_partial_window_call.val()
            + vta.reliability_derated_vta_npv.val() * 1e-2
            + fl.constrained_information_ratio.val()
            + acd.duration_adjusted_price_impact_bps.val() * 1e-2
            + hlz.deflated_sharpe_ratio_dsr.val()
            + ml.mlmc_telescoping_price.val()
            + inf
    }
    pub fn mf_rates_wildcard_ctd_pivot_tarf_bergomi_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        Self::mf_rates_wildcard_ctd_pivot_tarf_bergomi_live_metric_from_market(&market, bump)
    }
    pub fn kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let bump = if matches!(self.scenario, ScenarioKind::CreditWideningPlus100Bp) {
            0.01
        } else {
            0.0
        };
        Self::kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_live_metric_from_market(&market, bump)
    }
    pub fn mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_live_metric(&self) -> f64 {
        let market = self.scenario_market();
        let bump = if matches!(self.scenario, ScenarioKind::CreditWideningPlus100Bp) {
            0.01
        } else {
            0.0
        };
        Self::mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_live_metric_from_market(&market, bump)
    }
    pub fn step_tick(&mut self) -> WorkstationFrameData {
        self.tick += 1;
        let market = self.scenario_market();
        let start = Instant::now();
        let summary = self.dag.evaluate_portfolio(&market, &self.trades);
        let elapsed = start.elapsed().as_secs_f64() * 1e6;
        let (credit_raw, jtd, _, cds_npv) = self.credit_metrics(self.scenario);
        let (cva, dva, fva) = self.xva_metrics();
        let mva = self.mva_metric();
        let (inflation, inflation_npv) = self.inflation_metrics(self.scenario);
        let active_curves = self.active_curve_system(&market);
        let curve = FlatForwardCurve::new(self.base_market.ref_date, market.ois_quotes[0]);
        let credit_bump = if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
            0.01
        } else {
            0.0
        };
        // Live telemetry evaluates representative engines on the
        // current market/scenario snapshot so every tick and stress path is
        // reflected in the workstation metric.
        let sinhou_metric = price_sinh_ou_rfr_caplet(
            &SinhOuRfrCapletSpec {
                notional: 100.0,
                start_time: 0.5,
                end_time: 0.75,
                strike: market.ois_quotes[0] + credit_bump,
                alpha_mean_reversion: 0.05,
                sigma_short_rate: (0.01 + market.eq_vol * 0.03).max(1e-4),
                gamma_smile: 1.0 + market.eq_vol,
                y_star_skew: -0.15,
                is_caplet: true,
                is_backward_compounded: true,
            },
            &curve,
        )
        .npv
        .val();
        let mrb_metric = evaluate_mixed_rough_bergomi_vix(&MixedRoughBergomiVixSpec {
            xi0: (market.eq_vol * market.eq_vol).max(1e-4),
            hurst: 0.12,
            mixing_weight_lambda: 0.55,
            nu1: 0.9,
            nu2: 2.0 + market.eq_vol,
            expiry: 0.25,
            vix_window: 1.0 / 12.0,
            strike: 20.0,
            discount_factor: (-market.ois_quotes[0] * 0.25).exp(),
        })
        .vix_call_price
        .val();
        let choi_metric =
            ql_exotics::price_choi_sum_of_bsm_basket_option(&ql_exotics::ChoiSumOfBsmBasketSpec {
                spots: vec![market.eq_spot.max(1e-4), (0.9 * market.eq_spot).max(1e-4)],
                weights: vec![1.0, -1.0],
                vols: vec![market.eq_vol.max(1e-4), (0.8 * market.eq_vol).max(1e-4)],
                div_yields: vec![0.01, 0.02],
                correlation_matrix: vec![vec![1.0, 0.5], vec![0.5, 1.0]],
                strike: 5.0 + credit_bump * 100.0,
                risk_free_rate: market.ois_quotes[0],
                expiry: 1.0,
                is_call: true,
            })
            .price
            .val();
        let afs_metric = solve_alfonsi_fruth_schied_lob_execution(&AlfonsiFruthSchiedLobSpec {
            total_order_size: 100.0 + market.eq_spot * 0.01,
            horizon: 1.0,
            num_intervals: 4,
            lob_density_q: 50.0,
            shape_exponent_alpha: 0.4,
            resilience_rho: 1.2 + market.eq_vol,
        })
        .volume_resilience_cost
        .val();
        let exvar_metric =
            evaluate_and_optimize_expectile_portfolio(&ExpectileRiskAndPortfolioSpec {
                scenario_returns: vec![
                    vec![-0.07 - credit_bump, -0.01],
                    vec![-0.04, -0.005],
                    vec![0.02 - market.eq_vol * 0.01, 0.008],
                    vec![0.05, 0.012],
                ],
                initial_weights: vec![0.5, 0.5],
                tau_confidence: 0.85,
                return_weight_lambda: 0.2,
                smoothing_eps: 1e-3,
            })
            .optimal_expectile_exvar
            .val();
        let speed_metric =
            evaluate_third_order_speed_vol_expansion(&ThirdOrderSpeedVolExpansionSpec {
                xi0: (market.eq_vol * market.eq_vol).max(1e-4),
                omega_vol_of_vol: 1.8,
                k_mean_reversion: 3.0,
                expiry: 0.25,
                window_delta: 1.0 / 12.0,
                strike: 0.20,
                discount_factor: (-market.ois_quotes[0] * 0.25).exp(),
            })
            .call_price
            .val();

        // Complete telemetry: all fourteen engines consume the live
        // market/scenario snapshot so rates, equity, credit and tick changes
        // propagate through the aggregate metric.
        let cir_metric = ql_rates::price_cirpp_sofr_futures(
            &ql_rates::CirppSofrFuturesSpec {
                start_time: 0.5,
                end_time: 0.75,
                x0: market.ois_quotes[0] + credit_bump,
                b_drift: 0.006,
                beta_mean_reversion: 0.20,
                sigma: 0.06 + market.eq_vol * 0.02,
            },
            &curve,
        )
        .futures_rate_3m_compounded
        .val();
        let bk_metric = ql_fixed_income::analyze_black_karasinski_callable_frn_oas_ad(
            &ql_fixed_income::BlackKarasinskiCallableFrnSpec {
                face_amount: 100.0,
                quoted_margin: 0.008,
                cap_rate: 0.055,
                floor_rate: 0.020,
                maturity: 2.0,
                coupon_steps: 6,
                mean_reversion: 0.15,
                log_rate_vol: 0.22 + market.eq_vol * 0.005,
                flat_zero_rate: market.ois_quotes[0],
                call_price: 100.2,
                first_call_time: 1.0,
                target_market_price: 99.8 + credit_bump * 100.0,
                smoothing_eps: 1e-3,
            },
        )
        .oas_spread
        .val();
        let surv = ql_credit::SurvivalProbabilityCurve::new(
            self.base_market.ref_date,
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            vec![0.025 + credit_bump; 5],
        );
        let mb_metric = ql_credit::price_morini_brigo_credit_index_option(
            &ql_credit::MoriniBrigoIndexOptionSpec {
                notional: 100.0,
                strike_spread: 0.025 + credit_bump,
                option_expiry: 1.0,
                payment_times: vec![2.0, 3.0, 4.0, 5.0],
                recovery_rate: 0.40,
                spread_vol: 0.45 + market.eq_vol * 0.01,
                num_constituents: 50,
                systemic_copula_corr: 0.60,
                is_payer: true,
            },
            &curve,
            &surv,
        )
        .npv
        .val();
        let wwm_metric = ql_xva::compute_brigo_vrins_wwm_cva(&ql_xva::BrigoVrinsWwmCvaSpec {
            times: vec![0.5, 1.0, 2.0],
            discount_rate: market.ois_quotes[0],
            base_mtm_profile: vec![10.0 + market.eq_spot * 0.01, 12.0, 11.0],
            exposure_vol_nu: 8.0,
            is_brownian_bridge_swap: false,
            y0: 0.02,
            cir_kappa: 0.35,
            cir_theta: 0.025,
            cir_vol: 0.08,
            intensity_shift_psi: 0.005,
            wwr_correlation: 0.60,
            lgd: 0.60,
        })
        .wwm_cva
        .val();
        let bcs_metric =
            ql_xva::compute_bichuch_capponi_sturm_xva_band(&ql_xva::AsymmetricRateXvaBandSpec {
                maturity: 1.0,
                num_steps: 6,
                clean_claim_value: 25.0 + market.eq_spot * 0.01,
                clean_delta_notional: 40.0,
                asset_vol: 0.20 + market.eq_vol * 0.005,
                collateral_alpha: 0.50,
                dividend_discount_rate: 0.02,
                repo_lend_rate: 0.020,
                repo_borrow_rate: 0.024,
                funding_lend_rate: 0.022,
                funding_borrow_rate: 0.036 + credit_bump,
                collateral_lend_rate: 0.020,
                collateral_borrow_rate: 0.022,
                own_bond_yield_mu_i: 0.042,
                cpty_bond_yield_mu_c: 0.048 + credit_bump,
                own_lgd: 0.60,
                cpty_lgd: 0.60,
            })
            .no_arbitrage_band_width
            .val();
        let ldz_metric = ql_exotics::price_li_deng_zhou_multi_spread_option(
            &ql_exotics::LiDengZhouMultiSpreadSpec {
                long_spot: market.eq_spot,
                long_div_yield: 0.01,
                long_vol: market.eq_vol.max(1e-4),
                short_spots: vec![0.9 * market.eq_spot, 0.85 * market.eq_spot],
                short_weights: vec![0.5, 0.5],
                short_div_yields: vec![0.015, 0.02],
                short_vols: vec![0.20, 0.22],
                correlation_matrix: vec![
                    vec![1.0, 0.5, 0.4],
                    vec![0.5, 1.0, 0.6],
                    vec![0.4, 0.6, 1.0],
                ],
                strike: 5.0 + credit_bump * 100.0,
                risk_free_rate: market.ois_quotes[0],
                expiry: 1.0,
                is_call: true,
            },
        )
        .price
        .val();
        let mrjd_metric = ql_commodities::price_hikspoors_jaimungal_mrjd_commodity(
            &ql_commodities::HikspoorsJaimungalMrjdSpec {
                x0: 0.10,
                y0: 4.20 + market.eq_spot * 0.001,
                j0_spike: 0.04,
                seasonal_phi: 0.0,
                kappa_x: 1.5,
                sigma_x: 0.35 + market.eq_vol * 0.005,
                mu_y_rn: 0.02,
                sigma_y: 0.18,
                rho_xy: -0.2,
                kappa_j_spike: 12.0,
                jump_intensity: 2.5 + credit_bump * 10.0,
                mean_jump_size: 0.22,
                option_expiry: 0.5,
                futures_maturity: 1.0,
                strike: 70.0,
                risk_free_rate: market.ois_quotes[0],
                is_call: true,
            },
        )
        .futures_option_price
        .val();
        let pcrp_metric = ql_portfolio::solve_principal_component_risk_parity(
            &ql_portfolio::FactorEigenRiskParitySpec {
                covariance: vec![
                    vec![0.04 + market.eq_vol * 0.0005, 0.02],
                    vec![0.02, 0.05 + credit_bump * 0.01],
                ],
                factor_risk_budgets: vec![0.5, 0.5],
            },
        )
        .weights[0]
            .val();
        let sinhou_mb_rbergomi_choi_afs_expectile_metric = sinhou_metric
            + mrb_metric
            + choi_metric
            + afs_metric
            + exvar_metric
            + speed_metric
            + cir_metric
            + bk_metric
            + mb_metric
            + wwm_metric
            + bcs_metric
            + ldz_metric
            + mrjd_metric
            + pcrp_metric;
        let fx = FxVannaVolgaSmile::new(
            market.eq_spot,
            0.03,
            0.01,
            1.0,
            FxSmileQuote {
                atm_vol: 0.1,
                rr_25d: 0.01,
                bf_25d: 0.005,
            },
        );
        let fx_delta =
            garman_kohlhagen_spot_delta(market.eq_spot, fx.k3, fx.sigma3, 1.0, 0.03, 0.01, true)
                .val();
        // Assemble the actual cross-asset SIMM buckets from the current
        // sensitivities, rather than applying one arbitrary scalar weight to
        // the whole IR ladder.  All values are expressed in the same dollar
        // sensitivity units before the risk weights are applied.
        let mut buckets = vec![SimmBucket {
            sensitivities: summary
                .total_par_dv01
                .iter()
                .map(|x| ql_ad::TapeVar::from_f64(*x * 1e-4))
                .collect(),
            risk_weights: vec![50.0; summary.total_par_dv01.len()],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0; summary.total_par_dv01.len()]; summary.total_par_dv01.len()],
        }];
        buckets.push(SimmBucket {
            sensitivities: credit_raw
                .iter()
                .map(|(_, value)| ql_ad::TapeVar::from_f64(*value * 1e-4))
                .collect(),
            risk_weights: vec![100.0; credit_raw.len()],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0; credit_raw.len()]; credit_raw.len()],
        });
        buckets.push(SimmBucket {
            sensitivities: vec![ql_ad::TapeVar::from_f64(
                summary.total_eq_delta * market.eq_spot,
            )],
            risk_weights: vec![0.15],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0]],
        });
        buckets.push(SimmBucket {
            sensitivities: vec![ql_ad::TapeVar::from_f64(fx_delta * market.eq_spot)],
            risk_weights: vec![0.08],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0]],
        });
        buckets.push(SimmBucket {
            sensitivities: vec![ql_ad::TapeVar::from_f64(inflation)],
            risk_weights: vec![50.0],
            concentration_threshold: 1e6,
            intra_corr: vec![vec![1.0]],
        });
        let class = |bucket: &SimmBucket<ql_ad::TapeVar<f64>>| {
            compute_isda_simm(std::slice::from_ref(bucket), &[vec![1.0]])
        };
        let ir_margin = compute_isda_simm(
            &[buckets[0].clone(), buckets[4].clone()],
            &[vec![1.0, 0.5], vec![0.5, 1.0]],
        );
        let class_margins = [
            ir_margin,
            class(&buckets[1]),
            ql_ad::TapeVar::zero(),
            class(&buckets[2]),
            ql_ad::TapeVar::zero(),
            class(&buckets[3]),
        ];
        let simm = compute_isda_simm_6_risk_classes(&class_margins).val();
        let adjusted_total_npv = summary.total_npv + cds_npv + inflation_npv;
        let prev = self
            .npv_history
            .last()
            .copied()
            .unwrap_or(adjusted_total_npv);
        self.npv_history.push(adjusted_total_npv);
        self.last_summary = summary.clone();
        let rate = active_curves.curves[0].zero_rate(1.0).val();
        let eq_vega = {
            let params = HestonParams {
                v0: ql_ad::Dual::<1>::variable(market.eq_vol * market.eq_vol, 0),
                kappa: ql_ad::Dual::from_value(market.heston_params.kappa),
                theta: ql_ad::Dual::from_value(market.heston_params.theta),
                sigma_v: ql_ad::Dual::from_value(market.heston_params.sigma_v),
                rho: ql_ad::Dual::from_value(market.heston_params.rho),
            };
            let v = ql_models::heston_cos_price(
                ql_ad::Dual::from_value(market.eq_spot),
                ql_ad::Dual::from_value(100.0),
                ql_ad::Dual::from_value(1.0),
                ql_ad::Dual::from_value(rate),
                ql_ad::Dual::from_value(0.0),
                &params,
                128,
                10.0,
                true,
            );
            v.d[0] * 2.0 * market.eq_vol
        };
        let ir_buckets: Vec<(String, f64)> = summary
            .total_par_dv01
            .iter()
            .enumerate()
            .map(|(i, x)| {
                let label = ["OIS 1Y", "OIS 2Y", "OIS 5Y", "FWD 1Y", "FWD 2Y", "FWD 5Y"]
                    .get(i)
                    .map(|s| (*s).to_string())
                    .unwrap_or_else(|| format!("Pillar {}", i + 1));
                (label, *x * 1e-4)
            })
            .collect();
        let credit_buckets: Vec<(String, f64)> = credit_raw
            .iter()
            .map(|(label, value)| (label.clone(), *value * 1e-4))
            .collect();
        let irs_npv = if let (Some(ois), Some(fwd)) =
            (active_curves.curves.first(), active_curves.curves.get(1))
        {
            let spec = VanillaSwapSpec {
                effective_date: market.ref_date,
                termination_date: market.ref_date.add_months(60, false),
                notional: 1_000_000.0,
                fixed_rate_or_coupon: 0.02,
                fixed_tenor: Period::Years(1),
                float_tenor: Period::Years(1),
                calendar: Arc::new(BitsetCalendar::weekend_only()),
                bdc: BusinessDayConvention::Unadjusted,
                fixed_dc: DayCounter::Act365Fixed,
                float_dc: DayCounter::Act365Fixed,
                swap_type: SwapType::Payer,
            };
            spec.resolve(market.ref_date).npv(ois, fwd)
        } else {
            0.0
        };
        let heston_params = HestonParams {
            v0: market.heston_params.v0,
            kappa: market.heston_params.kappa,
            theta: market.heston_params.theta,
            sigma_v: market.heston_params.sigma_v,
            rho: market.heston_params.rho,
        };
        let heston_npv = ql_models::heston_cos_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &heston_params,
            128,
            10.0,
            true,
        );
        let script_npv = match self.trades.get(5) {
            Some(PortfolioTrade::ScriptedNote { program, .. }) => ql_scripting::ScriptMcEngine {
                n_paths: 1024,
                seed: 42,
            }
            .price_gbm(program, market.eq_spot, market.eq_vol, rate, 0.0)
            .val(),
            _ => 0.0,
        };
        let (_, inflation_curve) = self.inflation_calibration(self.scenario);
        let zciis = ZeroCouponInflationSwap {
            notional: 1_000_000.0,
            fixed_rate: 0.021,
            maturity_years: 5.0,
            target_month: Some(7),
            pay_fixed: true,
        };
        let zciis_npv = zciis
            .npv(&inflation_curve.discount, &inflation_curve.curve)
            .val();
        let g2_swaption = EuropeanSwaption {
            notional: 1_000_000.0,
            strike: 0.025,
            expiry: 1.0,
            swap_times: (1..=6).map(|x| x as f64).collect(),
            is_payer: true,
            settlement: SwaptionSettlement::Physical,
        };
        let g2 = G2PlusPlus {
            a: 0.1,
            sigma: 0.01,
            b: 0.3,
            eta: 0.015,
            rho: -0.75,
        };
        let g2_npv = g2
            .swaption_price(&active_curves.curves[0], &g2_swaption, 16)
            .val();
        let vg = VarianceGamma {
            sigma: 0.18,
            nu: 0.20,
            theta: -0.10,
        };
        let levy_npv = cos_levy_price(
            &vg,
            market.eq_spot,
            100.0,
            rate,
            0.0,
            1.0,
            ql_instruments::OptionType::Call,
            128,
        )
        .val();
        let conv_spec = ConvertibleBondSpec {
            spot: market.eq_spot,
            conversion_ratio: 1.0,
            face_value: 100.0,
            coupon_rate: 0.03,
            coupon_times: (1..=3).map(|x| x as f64).collect(),
            maturity: 3.0,
            rate,
            div_yield: 0.0,
            credit_spread: 0.02,
            vol: market.eq_vol,
            call_price: None,
            put_price: None,
        };
        let (conv_npv, _) = price_convertible_bond_tf(&conv_spec, 48);
        let rfr_spec = RfrCapletSpec {
            notional: 1_000_000.0,
            forward_rate: 0.045,
            strike: 0.04,
            tau_0: 0.5,
            tau_1: 1.0,
            accrual_fraction: 0.5,
            discount_tau_1: active_curves.curves[0].discount_factor(1.0).val(),
            is_caplet: true,
            rate_type: RfrRateType::BackwardLooking,
        };
        let rfr_npv =
            willems_backward_sabr_caplet_price(&rfr_spec, 0.08, 0.5, -0.35, 0.45, 1.0).val();
        let xccs = ConstNotionalXccs {
            dom_notional: 1_100_000.0,
            for_notional: 1_000_000.0,
            spot_fx: 1.1,
            payment_times: (1..=5).map(|x| x as f64).collect(),
            accrual_fractions: vec![1.0; 5],
            for_basis_spread: 0.001,
            dom_spread: 0.0,
            receive_foreign: true,
        };
        let xccs_npv = xccs
            .npv(
                &active_curves.curves[0],
                &active_curves.curves[1],
                &active_curves.curves[1],
            )
            .val();
        let es = gaussian_expected_shortfall(
            ql_ad::TapeVar::from_f64(summary.total_npv),
            ql_ad::TapeVar::from_f64(summary.total_eq_delta.abs().max(1.0)),
            0.975,
        )
        .val();
        let rough_params = LiftedHestonParams {
            v0: market.heston_params.v0,
            lambda: market.heston_params.kappa,
            theta: market.heston_params.theta,
            nu: market.heston_params.sigma_v,
            rho: market.heston_params.rho,
            hurst: 0.1,
            n_factors: 8,
            r_n: 2.5,
        };
        let rough_npv = lifted_heston_cos_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &rough_params,
            64,
            8.0,
            16,
            true,
        )
        .val();
        let barrier_spec = BarrierOptionSpec {
            option_type: ql_instruments::OptionType::Call,
            barrier_type: ExoticBarrierType::DownOut,
            spot: market.eq_spot,
            strike: 100.0,
            barrier: 90.0,
            rebate: 0.0,
            rate,
            div_yield: 0.0,
            vol: market.eq_vol,
            expiry: 1.0,
            monitoring: MonitoringFrequency::Discrete { observations: 12 },
        };
        let barrier_npv = price_barrier_option(&barrier_spec).val();
        let commodity = SchwartzSmith2F {
            chi0: 0.0,
            xi0: market.eq_spot.ln(),
            kappa: 1.0,
            sigma_chi: 0.3,
            mu_xi: 0.0,
            sigma_xi: 0.2,
            rho: -0.3,
        };
        let commodity_futures = commodity.futures_price(1.0);
        let commodity_option = commodity
            .price_futures_option(
                1.0,
                1.5,
                commodity_futures,
                rate,
                ql_instruments::OptionType::Call,
            )
            .val();
        let crack = kirk_spread_option(
            commodity_futures,
            3.0,
            0.0,
            0.3,
            0.25,
            0.2,
            rate,
            1.0,
            ql_instruments::OptionType::Call,
        )
        .val();
        // Live telemetry for the final two workstation tabs.  These
        // are deliberately computed from the active scenario curves rather
        // than being display-only constants.
        let bermudan = BermudanSwaption {
            notional: 1_000_000.0,
            strike: 0.02,
            exercise_times: vec![1.0, 2.0, 3.0],
            swap_times: (0..=5).map(|x| x as f64).collect(),
            is_payer: true,
            settlement: SwaptionSettlement::Physical,
        };
        let bermudan_value = HullWhiteBermudanTreeEngine { steps: 30 }
            .price_single_curve(
                &HullWhite1F {
                    a: 0.05,
                    sigma: 0.01,
                },
                &active_curves.curves[0],
                &bermudan,
            )
            .val();
        let sa_report = compute_sa_ccr_ead(&SaCcrNettingSet {
            trades: vec![SaCcrTrade {
                mtm: ql_ad::TapeVar::from_f64(irs_npv),
                notional: ql_ad::TapeVar::from_f64(1_000_000.0),
                maturity: 5.0,
                asset_class: SaCcrAssetClass::InterestRate {
                    ccy: "USD".to_string(),
                    start_time: 0.0,
                    end_time: 5.0,
                },
                is_long_primary: true,
                option_spec: None,
            }],
            collateral: ql_ad::TapeVar::zero(),
            is_margined: false,
            threshold: ql_ad::TapeVar::zero(),
            mta: ql_ad::TapeVar::zero(),
            nica: ql_ad::TapeVar::zero(),
            mpor_years: 0.5,
        });
        let commodity_curve = CommodityFuturesCurve::from_futures_pillars(
            vec![1.0, 2.0, 3.0],
            &[80.0, 81.0, 82.0],
            ql_commodities::SeasonalParams {
                a1: 0.0,
                b1: 0.0,
                a2: 0.0,
                b2: 0.0,
            },
        );
        let swing_value = price_swing_option(
            &SwingOptionSpec {
                exercise_times: vec![1.0, 2.0, 3.0],
                strike: 80.0,
                unit_volume: 1.0,
                min_exercises: 0,
                max_exercises: 2,
                is_call: true,
            },
            &commodity_curve,
            1.0,
            0.3,
            rate,
            2,
        )
        .val();
        let qp = box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            expected_returns: vec![0.06, market.ois_quotes[2]],
            risk_aversion: 3.0,
            lower_bounds: vec![0.15, 0.15],
            upper_bounds: vec![0.80, 0.80],
        });
        let qp_value = qp[0];
        let storage_value = price_gas_storage(
            &GasStorageSpec {
                decision_times: vec![1.0, 2.0],
                max_capacity: 10.0,
                initial_inventory: 5.0,
                max_injection_rate: 5.0,
                max_withdrawal_rate: 5.0,
                injection_loss: 0.0,
                withdrawal_loss: 0.0,
                injection_cost: 0.0,
                withdrawal_cost: 0.0,
                inventory_grid_points: 10,
            },
            &commodity_curve,
            1.0,
            0.3,
            rate,
        );
        let sabr_forward = market.eq_spot / 100.0;
        let sabr_quotes = (0..5)
            .map(|i| SmileStrikeQuote {
                strike: sabr_forward * (0.96 + 0.02 * i as f64),
                market_vol: 0.20 + 0.004 * (i as f64 - 2.0).abs(),
                weight: 1.0,
            })
            .collect::<Vec<_>>();
        let sabr_live =
            CalibratedSabrSmile::calibrate(sabr_forward, 1.0, 0.5, (0.1, -0.2, 0.4), &sabr_quotes)
                .expect("live SABR telemetry");
        let (sabr_atm, sabr_vega) = sabr_live.vol_and_quote_vega(sabr_forward);
        let sabr_metric = sabr_atm + sabr_vega.iter().sum::<f64>();
        let rough_metric = rbergomi_hybrid_mc_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &RoughBergomiParams {
                xi0: market.heston_params.v0,
                eta: market.heston_params.sigma_v,
                rho: market.heston_params.rho,
                hurst: 0.1,
            },
            8,
            64,
            true,
            self.tick as u64,
        )
        .val()
            + price_worst_of_autocallable(
                &WorstOfAutocallableSpec {
                    notional: 1.0,
                    spots: vec![market.eq_spot, market.eq_spot * 0.98],
                    vols: vec![market.eq_vol, market.eq_vol * 1.05],
                    divs: vec![0.0, 0.0],
                    rate,
                    correlation: vec![vec![1.0, 0.5], vec![0.5, 1.0]],
                    observation_times: vec![0.5, 1.0],
                    autocall_barriers: vec![1.05, 1.05],
                    redemption_bonuses: vec![0.0, 0.0],
                    coupon_barriers: vec![0.8, 0.8],
                    coupon_rates: vec![0.03, 0.03],
                    memory_coupon: true,
                    knock_in_barrier: 0.65,
                    put_strike: 1.0,
                    observe_ki_all_dates: true,
                    smoothing_half_width: 0.01,
                    smooth_min_temp: 0.01,
                },
                32,
                self.tick as u64 + 17,
            )
            .unwrap_or(0.0);
        let execution_metric = obizhaeva_wang_execution(&ObizhaevaWangSpec {
            initial_inventory: 1_000_000.0,
            trade_times: (0..=8).map(|i| i as f64 / 8.0).collect(),
            volatility: market.eq_vol,
            transient_impact: 0.01,
            resilience_speed: 0.5,
            temporary_impact: 0.001,
            risk_aversion: 0.1,
        })
        .mean_variance_utility
        .val()
            + compute_frtb_sa(&ql_xva::FrtbSaSpec {
                buckets: vec![ql_xva::FrtbSensitivityBucket {
                    risk_class: ql_xva::FrtbRiskClass::Girr,
                    bucket_id: 0,
                    delta_sensitivities: vec![ql_ad::TapeVar::from_f64(
                        summary.total_par_dv01.iter().sum::<f64>(),
                    )],
                    delta_risk_weights: vec![0.015],
                    vega_sensitivities: vec![ql_ad::TapeVar::from_f64(sabr_metric)],
                    vega_risk_weights: vec![0.55],
                    curvature_up_pnl: vec![ql_ad::TapeVar::zero()],
                    curvature_down_pnl: vec![ql_ad::TapeVar::zero()],
                    curvature_shock_weights: vec![0.015],
                    intra_bucket_corr: vec![vec![1.0]],
                }],
                inter_bucket_corr: vec![],
                drc_obligors: vec![],
                exotic_residual_notionals: vec![],
                other_residual_notionals: vec![],
            })
            .sbm_med
            .val()
                * 1e-3
            + sabr_metric * 1e-3;
        // Tab 23 telemetry must be sourced from the four engines named
        // by the tab, rather than reusing the SLV/FRTB aggregates above.  Keep
        // each contribution tied to the active scenario so a stress tick is
        // observable in the workstation frame and in the rendered dashboard.
        let monotone_curve = MonotoneConvexZeroCurve::new(
            market.ref_date,
            vec![1.0, 2.0, 5.0],
            market.ois_quotes.to_vec(),
        );
        let monotone_metric = [0.5, 1.5, 3.0]
            .iter()
            .map(|&t| monotone_curve.zero_rate(t).val())
            .sum::<f64>();
        let vix_params = BatesParams {
            heston: HestonParams {
                v0: market.heston_params.v0,
                kappa: market.heston_params.kappa,
                theta: market.heston_params.theta,
                sigma_v: market.heston_params.sigma_v,
                rho: market.heston_params.rho,
            },
            lambda: 0.1,
            mu_j: -0.05,
            sigma_j: 0.2,
        };
        let vix =
            heston_bates_vix_future_and_option_price(&vix_params, 1.0, 30.0 / 365.0, 0.2, 0.98);
        let parisian = ParisianBarrierOptionSpec {
            option_type: ql_instruments::OptionType::Call,
            barrier_type: ExoticBarrierType::DownOut,
            clock_kind: ExcursionClockKind::ParisianConsecutive,
            spot: market.eq_spot,
            strike: 100.0,
            barrier: 90.0,
            rebate: 0.0,
            window_duration: 0.20,
            rate,
            div_yield: 0.0,
            vol: market.eq_vol,
            expiry: 1.0,
        };
        let parisian_metric = price_parisian_option_mc_ad(&parisian, 64, 16, 0.01, 23).val();
        let robust_weights = robust_markowitz_weights(&RobustMarkowitzSpec {
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            return_uncertainty_cov: vec![vec![0.01, 0.0], vec![0.0, 0.02]],
            expected_returns: vec![0.06, market.ois_quotes[2]],
            risk_aversion: 3.0,
            uncertainty_radius: 0.1,
            lower_bounds: vec![0.0, 0.0],
            upper_bounds: vec![0.8, 0.8],
        });
        let robust_metric = robust_weights
            .iter()
            .enumerate()
            .map(|(i, weight)| (i + 1) as f64 * weight)
            .sum::<f64>();
        let vix_parisian_robust_metric = monotone_metric
            + vix.vix_future
            + vix.call_price
            + vix.put_price
            + parisian_metric
            + robust_metric;
        let dh_metric = double_heston_cos_price(
            market.eq_spot,
            100.0,
            1.0,
            rate,
            0.0,
            &DoubleHestonParams {
                factor1: heston_params,
                factor2: HestonParams {
                    v0: 0.01,
                    kappa: 2.0,
                    theta: 0.01,
                    sigma_v: 0.2,
                    rho: -0.3,
                },
            },
            48,
            10.0,
            true,
        )
        .val();
        let tarf_metric = price_fx_tarf_mc_ad(
            &FxTarfSpec {
                notional: 100.0,
                spot: 1.05,
                strike: 1.0,
                r_dom: rate,
                r_for: 0.01,
                vol: 0.12,
                fixing_times: vec![0.25, 0.5, 0.75, 1.0],
                target_profit: 0.15,
                downside_leverage: 1.5,
                is_exporter_call: true,
                cap_convention: FxTarfCapConvention::ExactCap,
                smoothing_eps: 1e-3,
            },
            32,
            7,
        )
        .npv;
        let coco_metric = price_coco_bond_ds(&CocoBondSpec {
            face_value: 100.0,
            coupon_rate: 0.05,
            coupon_times: vec![1.0, 2.0],
            maturity: 2.0,
            rate,
            credit_spread: 0.015,
            spot: market.eq_spot,
            trigger_barrier: 60.0,
            div_yield: 0.0,
            eq_vol: market.eq_vol,
            trigger_payoff: CocoTriggerPayoff::WriteDown {
                recovery_fraction: 0.4,
            },
        })
        .dirty_price;
        let hrp_w = hierarchical_risk_parity_weights(&[vec![0.04, 0.01], vec![0.01, 0.09]]);
        let double_heston_tarf_coco_hrp_metric =
            dh_metric + tarf_metric + coco_metric + hrp_w[0] + 2.0 * hrp_w[1];
        let jy_curve = ZeroInflationCurve::new(
            market.ref_date,
            100.0,
            vec![1.0, 2.0, 5.0, 10.0],
            vec![0.02, 0.021, 0.023, 0.025],
            None,
        );
        let jy_instrument = YoYInflationCapFloor {
            notional: 1.0,
            maturity_years: 10,
            cap_rate: Some(0.025),
            floor_rate: None,
            gearing: 1.0,
            spread: 0.0,
            cap_floor_type: YoYInflationCapFloorType::Cap,
        };
        let jy_value = price_yoy_capfloor_jarrow_yildirim(
            &jy_instrument,
            &active_curves.curves[0],
            &jy_curve,
            &JarrowYildirimParams {
                a_nom: 0.1,
                sigma_nom: 0.01,
                a_real: 0.08,
                sigma_real: 0.012,
                sigma_cpi: 0.01,
                rho_nom_real: 0.2,
                rho_nom_cpi: 0.1,
                rho_real_cpi: -0.2,
            },
        )
        .0
        .val();
        let tba = analyze_tba_dollar_roll(&TbaDollarRollSpec {
            front_clean_price: 99.0,
            back_clean_price: 98.5,
            wac: 0.04,
            pass_through_rate: 0.035,
            remaining_months: 120,
            front_accrued_fraction: 0.5,
            back_accrued_fraction: 0.5,
            roll_period_years: 1.0,
            assumed_cpr: 0.06,
            repo_rate: rate,
        });
        let pool = price_specified_pool_mbs(
            &SpecifiedPoolMbsSpec {
                original_balance: 100.0,
                wac: 0.04,
                pass_through_rate: 0.035,
                term_months: 120,
                seasoning_months: 30,
                cpr_min: 0.02,
                cpr_max: 0.10,
                refi_slope: 2.0,
                refi_inflection: 0.0,
                primary_secondary_spread: 0.01,
                cdr: 0.002,
                loss_severity: 0.4,
            },
            &active_curves.curves[0],
            0.0,
        );
        let bg = bergomi_2f_implied_vol_bg(
            &Bergomi2FParams {
                xi0: 0.04,
                omega: 0.4,
                theta: 0.5,
                k1: 1.0,
                k2: 2.0,
                rho12: 0.25,
                rho_s1: -0.5,
                rho_s2: -0.3,
            },
            100.0,
            100.0,
            1.0,
        )
        .val();
        let returns = (0..16)
            .map(|i| {
                vec![
                    0.01 * (i as f64 - 7.0),
                    0.02 * ((i % 5) as f64 - 2.0),
                    0.015 * ((i % 3) as f64 - 1.0),
                    0.01 * ((i % 4) as f64 - 1.5),
                ]
            })
            .collect::<Vec<_>>();
        let gerber = gerber_covariance_matrix(&returns, 0.5, 1e-3).unwrap_or_default();
        let nco =
            nested_clustered_optimization_weights(&ql_portfolio::NestedClusteredOptimizationSpec {
                covariance: gerber.clone(),
                expected_returns: None,
                risk_aversion: 1.0,
                num_clusters: 2,
            });
        let jy_tba_bergomi_nco_metric = jy_value
            + tba.net_roll_advantage.val()
            + pool.pv.val() * 1e-3
            + bg
            + nco.iter().sum::<f64>();
        let cheyette = price_cheyette_swaption_markov_projection(
            &CheyetteSwaptionSpec {
                notional: 100.0,
                strike: 0.032,
                option_expiry: 1.0,
                payment_times: vec![2.0, 3.0, 4.0, 5.0],
                mean_reversion: 0.05,
                sigma0: 0.008,
                skew_b: -1.5,
                is_payer: true,
            },
            &active_curves.curves[0],
        );
        let cmo = price_cmo_waterfall(
            &CmoWaterfallSpec {
                original_balance: 100.0,
                wac: 0.05,
                pass_through_rate: 0.045,
                term_months: 24,
                psa_speed: 150.0,
                pac_collar_psa: Some((100.0, 300.0)),
                tranches: vec![
                    CmoTrancheSpec {
                        initial_balance: 60.0,
                        coupon_rate: 0.04,
                        is_pac: true,
                    },
                    CmoTrancheSpec {
                        initial_balance: 40.0,
                        coupon_rate: 0.045,
                        is_pac: false,
                    },
                ],
            },
            &active_curves.curves[0],
            0.0,
        );
        let qrh = qrh_implied_vol_expansion(
            &QuadraticRoughHestonParams {
                hurst: 0.1,
                nu: 0.35,
                lambda: 1.2,
                a: 0.45,
                b: 0.2,
                c: market.eq_vol * market.eq_vol,
                z0: 0.0,
            },
            market.eq_spot,
            market.eq_spot,
            1.0,
        );
        let ctd = analyze_treasury_futures_ctd(&TreasuryFuturesBasketSpec {
            futures_price: 98.5,
            delivery_time: 0.25,
            repo_rate: rate,
            notional_standard_yield: 0.06,
            deliverables: vec![
                TreasuryFutureDeliverable {
                    clean_price: 99.2,
                    coupon_rate: 0.06,
                    remaining_semiannual_periods: 16,
                    accrued_interest_spot: 0.5,
                    accrued_interest_delivery: 1.0,
                    interim_coupon: 0.0,
                    yield_vol: 0.008,
                    modified_duration: 6.5,
                },
                TreasuryFutureDeliverable {
                    clean_price: 89.0,
                    coupon_rate: 0.045,
                    remaining_semiannual_periods: 18,
                    accrued_interest_spot: 0.4,
                    accrued_interest_delivery: 0.8,
                    interim_coupon: 0.0,
                    yield_vol: 0.008,
                    modified_duration: 7.2,
                },
            ],
        })
        .option_adjusted_futures_price;
        let prdc = price_prdc_note_hybrid(&PrdcNoteSpec {
            notional: 100.0,
            spot_fx: 150.0,
            r_dom: 0.01,
            r_for: 0.04,
            a_dom: 0.05,
            sigma_dom: 0.008,
            a_for: 0.08,
            sigma_for: 0.01,
            sigma_fx: 0.12,
            rho_fx_dom: 0.3,
            rho_fx_for: -0.2,
            rho_dom_for: 0.25,
            coupon_times: vec![1.0, 2.0],
            foreign_gearing: 0.15,
            domestic_gearing: 0.10,
            coupon_cap: 0.08,
            fx_knockout_barrier: Some(180.0),
        })
        .npv;
        let factor_attr = factor_risk_attribution(
            &[0.4, 0.6],
            &[vec![1.1, 0.3], vec![0.8, 0.9]],
            &[vec![0.04, 0.005], vec![0.005, 0.02]],
            &[0.01, 0.015],
        )
        .total_volatility;
        let dd_returns = (0..16)
            .map(|i| -0.045 + 0.006 * i as f64 - 0.002 * rate)
            .collect::<Vec<_>>();
        let evar_cdar = conditional_drawdown_at_risk(&dd_returns, 0.75).cdar
            + entropic_value_at_risk(&dd_returns, 0.90);
        let base_quotes = vec![
            BaseCorrelationTrancheQuote {
                attachment: 0.0,
                detachment: 0.03,
                running_spread: 0.015 + rate,
                upfront: 0.0,
            },
            BaseCorrelationTrancheQuote {
                attachment: 0.03,
                detachment: 0.07,
                running_spread: 0.02 + rate,
                upfront: 0.0,
            },
        ];
        let base_corr = bootstrap_base_correlation_curve(
            &base_quotes,
            0.02 + rate * 0.1,
            0.4,
            rate,
            &[1.0, 2.0, 3.0, 4.0, 5.0],
        );
        let (_, base_fair, _) = ql_credit::price_mezzanine_tranche_base_correlation(
            &base_corr,
            0.02,
            0.05,
            0.02 + rate,
            100.0,
        );
        let cheyette_cmo_qrh_evar_metric = cheyette.npv.val()
            + cmo.total_structure_pv.val()
            + qrh.val()
            + ctd
            + prdc
            + evar_cdar
            + factor_attr
            + base_fair.val();
        let suite29_curve = FlatForwardCurve::new(market.ref_date, rate);
        let asw = analyze_bond_asset_swap_and_z_spread(
            &BondAssetSwapSpec {
                clean_price: 0.97,
                accrued_interest: 0.01,
                coupon_rate: 0.04,
                payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
                accrual_fractions: vec![1.0; 5],
            },
            &suite29_curve,
        );
        let mf = price_markov_functional_1f(
            &MarkovFunctionalSwaptionSpec {
                expiry: 1.0,
                payment_times: vec![2.0, 3.0, 4.0, 5.0],
                mean_reversion: 0.05 + rate * 0.1,
                atm_normal_vol: 0.008 + market.eq_vol * 0.01,
                sabr_rho: -0.25,
                sabr_nu: 0.35,
                strike: 0.031 + rate * 0.05,
                notional: 100.0,
                is_payer: true,
            },
            &suite29_curve,
            12,
        );
        let ba_cva = compute_ba_cva_capital(&BaCvaPortfolioSpec {
            counterparties: vec![BaCvaCounterpartySpec {
                risk_weight: 0.05,
                effective_maturity: 3.0,
                ead: 100.0 + market.eq_vol * 20.0,
                single_name_hedges: vec![BaCvaSingleNameHedge {
                    notional: 35.0,
                    maturity: 3.0,
                    hedge_correlation: 0.8,
                }],
            }],
            index_hedges: vec![],
            hurdle_rate: 0.10,
            risk_free_rate: rate,
        });
        let black_cox = analyze_structural_credit_black_cox(&BlackCoxStructuralSpec {
            asset_value: 140.0,
            debt_face_value: 100.0,
            default_barrier: 85.0,
            barrier_growth_rate: 0.02,
            asset_vol: 0.20 + market.eq_vol,
            risk_free_rate: rate + credit_bump,
            asset_payout_rate: 0.01,
            maturity: 3.0,
            recovery_rate: 0.40,
        });
        let ah_slice = solve_andreasen_huge_step(
            &AndreasenHugeGridSpec {
                forward: 100.0,
                rate,
                dt: 0.5,
                strike_min: 70.0,
                strike_max: 130.0,
                n_grid: 25,
                local_vol_atm: 0.18 + market.eq_vol,
                local_vol_skew: -0.15,
                local_vol_smile: 0.10,
            },
            None,
        );
        let ah_call = interpolate_andreasen_huge_call(&ah_slice, 100.0);
        let suite29_curve_futures = CommodityFuturesCurve::from_futures_pillars(
            vec![0.25, 0.50, 0.75, 1.0],
            &[
                45.0 + market.eq_vol * 5.0 + market.eq_spot * 1e-3,
                58.0,
                42.0,
                62.0,
            ],
            ql_commodities::SeasonalParams {
                a1: 0.0,
                b1: 0.0,
                a2: 0.0,
                b2: 0.0,
            },
        );
        let bess = price_bess_arbitrage_dp(
            &BessSpec {
                decision_times: vec![0.25, 0.50, 0.75, 1.0],
                energy_capacity_mwh: 4.0,
                initial_soc_mwh: 2.0,
                max_charge_mw: 8.0,
                max_discharge_mw: 8.0,
                charge_efficiency: 0.92,
                discharge_efficiency: 0.92,
                cycle_degradation_cost: 1.5,
                ancillary_reserve_rate: 2.0,
                soc_grid_points: 3,
            },
            &suite29_curve_futures,
            1.5,
            0.35 + market.eq_vol,
            rate,
        );
        let fx_accumulator = price_fx_accumulator_strip_no_ko(&FxAccumulatorSpec {
            base_notional_per_fixing: 1.0,
            spot: market.eq_spot,
            strike: 100.0,
            knockout_barrier: 120.0,
            r_dom: rate + credit_bump,
            r_for: 0.01,
            vol: market.eq_vol.max(0.01),
            fixing_times: vec![0.25, 0.50, 0.75, 1.0],
            guaranteed_fixings: 0,
            downside_gearing: 1.2,
            is_accumulator: true,
            smoothing_eps: 1e-6,
        });
        let atlas = price_atlas_option_mc_ad(
            &AtlasOptionSpec {
                spots: vec![market.eq_spot, market.eq_spot * 1.02, market.eq_spot * 0.98],
                vols: vec![market.eq_vol.max(0.01); 3],
                div_yields: vec![0.01; 3],
                correlation: 0.35,
                rate: rate + credit_bump,
                expiry: 1.0,
                strike: 0.0,
                notional: 1.0,
                remove_best: 0,
                remove_worst: 0,
            },
            32,
            17,
        );
        let take_or_pay = price_take_or_pay_contract(
            &TakeOrPayContractSpec {
                decision_times: vec![0.25, 0.50, 0.75, 1.0],
                strike: 50.0 + credit_bump,
                daily_min_qty: 1.0,
                daily_max_qty: 3.0,
                cumulative_min_qty: 6.0,
                cumulative_max_qty: 10.0,
                shortfall_penalty_rate: 50.0,
                daily_discretization_steps: 2,
            },
            &suite29_curve_futures,
            1.5,
            0.30 + market.eq_vol,
            rate + credit_bump,
        );
        let bl = idzorek_black_litterman_posterior(&IdzorekBlackLittermanSpec {
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            market_weights: vec![0.6, 0.4],
            risk_aversion: 2.5,
            tau: 0.05,
            view_matrix: vec![vec![1.0, -1.0]],
            view_returns: vec![0.03 + rate * 0.1],
            view_confidences: vec![0.65 + market.eq_vol * 0.1],
        });
        let kelly = kelly_fractional_allocation(&KellyAllocationSpec {
            expected_returns: vec![0.05 + rate, 0.07 + market.eq_vol],
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            risk_free_rate: rate,
            max_kelly_fraction: 0.50,
            max_drawdown_limit: 0.20,
            drawdown_tail_prob: 0.10,
        });
        let dark_lit = optimal_dark_lit_execution(&DarkLitRoutingSpec {
            initial_inventory: 100.0,
            horizon: 1.0,
            num_steps: 4,
            volatility: market.eq_vol,
            lit_temporary_impact: 1e-3,
            lit_permanent_impact: 2e-4,
            dark_fill_intensity: 1.5 + rate * 2.0,
            risk_aversion: 0.25,
        });
        let asw_mf1f_ah_bess_ep_metric = asw.z_spread.val()
            + asw.par_asset_swap_spread.val()
            + mf.swaption_npv.val()
            + mf.cms_swaplet_rate.val()
            + ba_cva.capital_charge.val()
            + black_cox.black_cox_cds_spread.val()
            + ah_call.val()
            + bess.total_npv.val()
            + bl.weights.iter().map(|x| x.val()).sum::<f64>()
            + kelly.effective_kelly_fraction.val()
            + kelly.expected_log_growth_rate.val()
            + dark_lit.expected_total_cost.val()
            + fx_accumulator.val()
            + atlas.val()
            + take_or_pay.npv.val()
            + market.eq_spot * 1e-6;
        let suite30_for_curve =
            FlatForwardCurve::new(market.ref_date, rate * 0.85 + inflation * 0.2);
        let suite30_csa_curve =
            FlatForwardCurve::new(market.ref_date, rate * 0.90 + credit_bump * 0.3);
        let mtm_xccs = price_mtm_resetting_xccs(
            &MtmResettingXccs {
                for_notional: 100.0,
                spot_fx: market.eq_spot + (market.eq_spot - 100.0) * 0.001,
                payment_times: vec![0.5, 1.0, 1.5, 2.0],
                accrual_fractions: vec![0.5; 4],
                for_basis_spread: -0.0015 - credit_bump * 0.2,
                dom_spread: 0.0,
                receive_foreign: true,
            },
            &suite29_curve,
            &suite30_for_curve,
            &suite30_csa_curve,
        );
        let tenor_basis = price_tenor_basis_swap(
            &TenorBasisSwapSpec {
                notional: 100.0,
                short_leg_times: vec![0.5, 1.0],
                long_leg_times: vec![1.0],
                short_leg_basis_spread: 0.001 + credit_bump * 0.1,
                pay_short_leg: true,
            },
            &suite29_curve,
            &suite30_for_curve,
            &suite30_csa_curve,
        );
        let callable_oas = analyze_callable_bond_oas_ad(&CallableBondOasSpec {
            market_dirty_price: 100.0 - credit_bump * 50.0,
            face_amount: 100.0,
            coupon_rate: 0.04 + inflation * 0.2,
            payment_times: vec![1.0, 2.0, 3.0],
            call_price: Some(101.0),
            put_price: None,
            short_rate: rate,
            short_rate_vol: 0.012 + market.eq_vol * 0.02,
            tree_steps: 8,
            smoothing_eps: 1e-3,
        });
        let ctd_colva = compute_multicurrency_ctd_colva(&MultiCurrencyCsaSpec {
            eligible_currency_spreads: vec![
                0.008 + rate * 0.1,
                0.010 + credit_bump * 0.5,
                0.009 + inflation * 0.3,
            ],
            spread_vols: vec![0.006 + market.eq_vol * 0.01, 0.008, 0.007],
            cross_currency_corr: 0.25,
            collateral_balance_profile: vec![50.0, 40.0],
            payment_times: vec![0.5, 1.0],
        });
        let turlakov_wwr = compute_turlakov_wwr_cva_fva(&TurlakovWwrSpec {
            times: vec![0.5, 1.0],
            base_epe: vec![10.0, 14.0],
            stressed_epe: vec![25.0 + (100.0 - market.eq_spot).max(0.0), 32.0],
            base_mtm: vec![8.0, 10.0],
            stressed_mtm: vec![20.0, 24.0],
            sovereign_hazard_rate: 0.02 + credit_bump,
            counterparty_hazard_rate: 0.025 + credit_bump,
            crisis_prob_threshold: 0.08,
            coupling_lambda: 0.45,
            lgd: 0.60,
            base_funding_spread: 0.01,
            stressed_funding_spread: 0.025 + credit_bump,
            discount_rate: rate,
        });
        let zabr_q = evaluate_zabr_smile(
            &ZabrParams {
                alpha: 0.18 + market.eq_vol * 0.5,
                beta: 0.7,
                nu: 0.40,
                rho: -0.25,
                gamma: 0.85,
            },
            market.eq_spot,
            105.0,
            1.0,
        );
        let pdv_iv = pdv_implied_vol_expansion(
            &GuyonPdvParams {
                beta0: 0.15 + market.eq_vol * 0.5,
                beta1: -0.12,
                beta2: 0.25,
                theta1: 0.5,
                theta2: 0.5,
                r10: -0.05,
                r11: -0.02,
                r20: 0.04,
                r21: 0.04,
            },
            market.eq_spot,
            95.0,
            1.0,
        );
        let fx_kiko = price_fx_kiko_option(&FxKikoOptionSpec {
            spot: market.eq_spot,
            strike: 100.0,
            lower_barrier: 75.0,
            upper_barrier: 125.0,
            vol: market.eq_vol,
            expiry: 1.0,
            r_dom: rate,
            r_for: 0.01,
            is_call: true,
            lower_is_ko_upper_is_ki: true,
        });
        let fx_wb = price_fx_window_barrier_option(&FxWindowBarrierSpec {
            spot: market.eq_spot,
            strike: 100.0,
            barrier: 130.0,
            vol: market.eq_vol,
            expiry: 1.0,
            r_dom: rate,
            r_for: 0.01,
            window_start: 0.1,
            window_end: 0.8,
            is_up_and_out: true,
            is_call: true,
        });
        let everest_val = price_everest_option_mc_ad(
            &EverestOptionSpec {
                spots: vec![market.eq_spot, 102.0, 98.0],
                vols: vec![market.eq_vol, 0.22, 0.20],
                correlation: 0.4,
                rate,
                expiry: 1.0,
                strike: 0.85,
                notional: 10.0,
                tau_temp: 0.02,
            },
            16,
        );
        let capped_timer = price_capped_timer_option_quadrature(
            &FiniteCapTimerOptionSpec {
                spot: market.eq_spot,
                strike: 100.0,
                rate,
                div_yield: 0.01,
                variance_budget: 0.04,
                t_max: 0.8,
                is_call: true,
            },
            &market.heston_params,
        );
        let weather_opt = price_weather_hdd_cdd_option(&WeatherDerivativeSpec {
            kind: DegreeDayKind::Hdd,
            temperature_base: 18.0,
            initial_temperature: 12.0 + inflation * 50.0,
            mean_level: 13.0,
            seasonal_amplitude: 5.0,
            seasonal_phase: 0.2,
            mean_reversion: 0.25,
            daily_vol: 2.0,
            market_price_risk: 0.05,
            times: vec![1.0, 2.0, 3.0],
            strike: 14.0,
            notional: 5.0,
            is_call: true,
            cap: Some(60.0),
        });
        let clean_spread = price_clean_spark_dark_spread(&CleanSparkDarkSpreadSpec {
            power: 80.0 + inflation * 100.0,
            fuel: 3.0,
            carbon: 20.0,
            heat_rate: 7.0,
            emission_factor: 0.5,
            strike: 10.0,
            vol_power: 0.22 + market.eq_vol * 0.2,
            vol_fuel: 0.30,
            vol_carbon: 0.35,
            rho_power_fuel: 0.45,
            rho_power_carbon: 0.25,
            rho_fuel_carbon: 0.20,
            rate,
            expiry: 1.0,
            option_type: OptionType::Call,
        });
        let min_torsion_enb = effective_number_of_bets_minimum_torsion(
            &[0.6, 0.4],
            &[
                vec![0.04 + market.eq_vol * 0.05, 0.012],
                vec![0.012, 0.09 + credit_bump * 0.5],
            ],
        );
        let ep_post = entropy_pooling_posterior(&EntropyPoolingSpec {
            scenario_returns: vec![vec![0.08, 0.02], vec![-0.04, 0.06], vec![0.02, -0.03]],
            prior_probs: vec![1.0 / 3.0; 3],
            view_matrix: vec![vec![0.08, -0.04, 0.02]],
            view_targets: vec![0.03 + rate * 0.2],
            risk_aversion: 2.5,
        });
        let mtmxccs_colva_zabr_weather_enb_metric = mtm_xccs.npv
            + tenor_basis.npv
            + callable_oas.oas * 100.0
            + ctd_colva.total_ctd_colva
            + turlakov_wwr.wwr_cva
            + zabr_q.call_price
            + pdv_iv
            + fx_kiko
            + fx_wb
            + everest_val
            + capped_timer.npv
            + weather_opt.npv
            + clean_spread
            + min_torsion_enb.enb
            + ep_post.weights[0];
        // Institutional engines: evaluate each live path against the
        // current scenario market and aggregate into the workstation telemetry.
        let suite31_mbs = analyze_mbs_burnout_oas_ad(
            &MbsBurnoutOasSpec {
                market_price: 101.0 + market.eq_spot * 0.01,
                original_balance: 100.0,
                wac: 0.06,
                pass_through_rate: 0.055,
                term_months: 24,
                seasoning_months: 12,
                turnover_cpr: 0.06,
                max_refi_cpr: 0.30,
                refi_slope: 60.0,
                refi_incentive_threshold: 0.005,
                burnout_beta: 10.0 + market.eq_vol,
                primary_secondary_spread: 0.01,
            },
            &active_curves.curves[0],
        )
        .oas;
        let suite31_zciis = price_zciis_capfloor_jarrow_yildirim(
            &ZciisCapFloorSpec {
                notional: 100.0,
                maturity_years: 5.0,
                strike_annual_rate: 0.02 + credit_bump,
                target_month: None,
                is_cap: true,
            },
            &active_curves.curves[0],
            &inflation_curve.curve,
            &JarrowYildirimParams {
                a_nom: 0.05,
                sigma_nom: 0.008,
                a_real: 0.04,
                sigma_real: 0.006,
                sigma_cpi: 0.012 + market.eq_vol * 0.01,
                rho_nom_real: 0.3,
                rho_nom_cpi: 0.2,
                rho_real_cpi: -0.1,
            },
        )
        .npv;
        let p19_credit = SurvivalProbabilityCurve::new(
            market.ref_date,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            vec![0.02 + credit_bump; 6],
        );
        let suite31_cmcds = price_cmcds_market_model(
            &CmCdsMarketModelSpec {
                notional: 100.0,
                grid_times: vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                start_idx: 0,
                end_idx: 3,
                cm_tenor_periods: 2,
                participation_rate: 0.95,
                lgd: 0.6,
                forward_cds_vols: vec![0.30 + market.eq_vol * 0.1; 6],
                forward_cds_corr: 0.8,
            },
            &active_curves.curves[0],
            &p19_credit,
        )
        .seller_npv;
        let suite31_simm = compute_dynamic_simm_mva_and_spike(&DynamicSimmMvaSpec {
            times: vec![0.5, 1.0, 1.5],
            bucket_sensitivities_t0: vec![100.0 + summary.total_eq_delta.abs(), 70.0],
            bucket_decay_rates: vec![0.25, 0.30],
            risk_weights: vec![0.5, 0.6],
            concentration_thresholds: vec![80.0, 80.0],
            inter_bucket_corr: 0.35,
            scheduled_cashflows: vec![12.0, 15.0, 18.0],
            mpor_years: 10.0 / 252.0,
            margin_funding_spread: 0.015 + credit_bump,
            hazard_rate: 0.02 + credit_bump,
            lgd: 0.6,
            discount_rate: rate,
        })
        .total_aps_cva;
        let suite31_jlt = evaluate_jlt_rating_migration_credit(&JltRatingMigrationSpec {
            generator_matrix: vec![
                vec![-0.06, 0.045, 0.01, 0.005],
                vec![0.02, -0.09, 0.055, 0.015],
                vec![0.005, 0.025, -0.11, 0.08],
                vec![0.0, 0.0, 0.0, 0.0],
            ],
            risk_premium_multipliers: vec![1.2 + credit_bump, 1.4, 1.6, 1.0],
            initial_rating_index: 0,
            horizon_years: 3.0,
            payment_times: vec![1.0, 2.0, 3.0],
            base_coupon_rate: 0.04,
            downgrade_step_up_spreads: vec![0.0, 0.005, 0.015, 0.0],
            recovery_rate: 0.4,
            risk_free_rate: rate,
        })
        .step_up_bond_price;
        let suite31_qou = evaluate_quintic_ou_spx_vix(
            &QuinticOuParams {
                xi0: 0.04 + market.eq_vol * 0.01,
                kappa: 2.5,
                nu: 0.8,
                rho: -0.7,
                alpha0: 1.0,
                alpha1: 0.30,
                alpha3: 0.10,
                alpha5: 0.02,
            },
            0.5,
            1.0 / 12.0,
            20.0,
            (-rate * 0.5).exp(),
        )
        .vix_call_price;
        let suite31_fbs = evaluate_free_boundary_sabr(
            &FreeBoundarySabrParams {
                forward: -0.003 + market.eq_spot * 1e-5,
                expiry: 1.0,
                alpha: 0.008 + market.eq_vol * 0.001,
                beta: 0.25,
                rho: -0.20,
                nu: 0.40,
            },
            -0.001,
            (-rate).exp(),
        )
        .call_price;
        let suite31_afv = price_convertible_bond_afv_pde(
            &AfvConvertibleBondSpec {
                spot: market.eq_spot,
                conversion_ratio: 1.0,
                face_value: 100.0,
                coupon_rate: 0.03,
                maturity: 1.5,
                rate,
                div_yield: 0.01,
                vol: market.eq_vol,
                base_hazard_rate: 0.02 + credit_bump,
                hazard_equity_elasticity: 1.2,
                equity_loss_on_default_eta: 0.5,
                bond_recovery_rate: 0.4,
                call_price: Some(115.0),
                soft_call_trigger_spot: Some(130.0),
                put_price: Some(95.0),
            },
            15,
            8,
        )
        .total_price;
        let suite31_gswap = price_gamma_swap_and_variance_option(&GammaSwapAndVarOptionSpec {
            spot: market.eq_spot,
            rate,
            div_yield: 0.01,
            expiry: 1.0,
            cut_strike: 100.0,
            variance_strike: 0.04,
            quotes: vec![
                VarianceSwapStrikeQuote {
                    strike: 85.0,
                    implied_vol: 0.24,
                },
                VarianceSwapStrikeQuote {
                    strike: 100.0,
                    implied_vol: 0.20,
                },
                VarianceSwapStrikeQuote {
                    strike: 115.0,
                    implied_vol: 0.18,
                },
            ],
            heston: HestonParams {
                v0: market.heston_params.v0,
                kappa: 2.0,
                theta: 0.04,
                sigma_v: 0.35,
                rho: -0.7,
            },
            is_variance_call: true,
        })
        .variance_option_price;
        let suite31_cs = price_clewlow_strickland_strip_swaption(&ClewlowStricklandSwaptionSpec {
            option_expiry: 0.5,
            delivery_times: vec![0.75, 1.0, 1.25],
            forward_prices: vec![75.0 + market.eq_spot * 0.01, 76.0, 74.0],
            delivery_weights: vec![1.0 / 3.0; 3],
            strike: 74.5,
            rate,
            sigma_short: 0.32 + market.eq_vol * 0.1,
            alpha_short: 1.1,
            sigma_long: 0.15,
            rho_factors: 0.2,
            seasonal_vol_amp: 0.05,
            is_payer_call: true,
        })
        .npv;
        let suite31_mm = solve_avellaneda_stoikov_market_making(&AvellanedaStoikovMmSpec {
            mid_price: market.eq_spot,
            initial_inventory: 3.0 + credit_bump * 100.0,
            horizon: 1.0,
            num_steps: 6,
            volatility: market.eq_vol,
            risk_aversion: 0.5,
            order_book_kappa: 1.5,
            base_arrival_intensity: 10.0,
            alpha_drift: 0.01,
        })
        .expected_spread_capture_pnl;
        let suite31_cdar = optimize_mean_cdar_portfolio(&MeanCdarPortfolioSpec {
            path_asset_returns: vec![
                vec![0.03, 0.01],
                vec![-0.06, -0.01],
                vec![-0.04, -0.005],
                vec![0.04, 0.012],
            ],
            expected_returns: vec![0.09 + market.eq_vol * 0.01, 0.04],
            cdar_aversion: 4.0,
            confidence: 0.80,
            smoothing_eps: 1e-3,
        })
        .cdar;
        let suite31_schur = schur_complementary_hrp(&SchurHrpSpec {
            covariance: vec![
                vec![0.04 + market.eq_vol * 0.01, 0.012],
                vec![0.012, 0.09 + credit_bump],
            ],
            gamma: 0.6,
        })
        .portfolio_variance;
        let suite31_oas = oas_and_nonlinear_shrinkage_covariance(&[
            vec![0.02 + market.eq_vol * 0.001, -0.01],
            vec![-0.015, 0.025 + credit_bump * 0.1],
            vec![0.01, 0.005],
            vec![-0.005, -0.01],
        ])
        .oas_shrinkage_intensity;
        let burnout_cmcds_quintic_afv_cs_schur_metric = suite31_mbs * 100.0
            + suite31_zciis
            + suite31_cmcds
            + suite31_simm
            + suite31_jlt
            + suite31_qou
            + suite31_fbs * 100.0
            + suite31_afv
            + suite31_gswap
            + suite31_cs
            + suite31_mm
            + suite31_cdar
            + suite31_schur
            + suite31_oas;
        let krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric =
            Self::krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_live_metric(&market, credit_bump);
        let blotter_rows = vec![
            ("IRS 5Y PayFix".to_string(), irs_npv, "IFT".to_string()),
            ("Heston Call 1Y".to_string(), heston_npv, "COS".to_string()),
            (
                "RoughHest H=.1".to_string(),
                rough_npv,
                "Lift 8-Fac".to_string(),
            ),
            ("Autocall Note".to_string(), script_npv, "AAD".to_string()),
            ("CDS 5Y BuyProt".to_string(), cds_npv, "IFT".to_string()),
            (
                "ZCIIS 5Y CPI".to_string(),
                zciis_npv,
                "Seasonal".to_string(),
            ),
            (
                "FX 25D VV Call".to_string(),
                fx.price_vanna_volga(fx.k3, true).val(),
                "Dual3".to_string(),
            ),
            (
                "EXOT-BAR-DN-OUT-1Y".to_string(),
                barrier_npv,
                "Dual5 BGK".to_string(),
            ),
            (
                "CMDTY-WTI-SS2F-1Y".to_string(),
                commodity_option + crack,
                "SS2F/Kirk".to_string(),
            ),
            ("G2++ Swaption".to_string(), g2_npv, "G2 Quad".to_string()),
            ("LEVY-VG-COS-1Y".to_string(), levy_npv, "COS VG".to_string()),
            (
                "CONV-TF-BOND-3Y".to_string(),
                conv_npv.val(),
                "TF + ASB".to_string(),
            ),
            (
                "SOFR RFR SABR 1Y".to_string(),
                rfr_npv,
                "Willems".to_string(),
            ),
            (
                "EURUSD XCCS 5Y".to_string(),
                xccs_npv,
                "IFT Basis".to_string(),
            ),
            ("Portfolio ES 97.5%".to_string(), es, "Acerbi".to_string()),
            (
                "Bermudan HW 3Y".to_string(),
                bermudan_value,
                "Tree IFT".to_string(),
            ),
            (
                "SA-CCR EAD".to_string(),
                sa_report.ead.val(),
                "BCBS279".to_string(),
            ),
            (
                "Commodity Swing".to_string(),
                swing_value,
                "OU DP".to_string(),
            ),
            ("Box Markowitz".to_string(), qp_value, "KKT IFT".to_string()),
            (
                "Gas Storage DP".to_string(),
                storage_value,
                "OU DP".to_string(),
            ),
        ];
        WorkstationFrameData {
            tick: self.tick,
            scenario: self.scenario,
            chart_mode: self.chart_mode,
            eval_micros: elapsed,
            dag_recomputed_trades: summary.last_trade_recomputes,
            dag_total_trades: self.trades.len(),
            dag_curves_rebuilt: summary.last_curve_recomputes > 0,
            total_npv: adjusted_total_npv,
            intraday_pnl: adjusted_total_npv - prev,
            eq_delta: summary.total_eq_delta,
            eq_vega,
            fx_spot_delta: fx_delta,
            inflation_par_delta_sum: inflation,
            cva,
            dva,
            fva,
            mva,
            simm_initial_margin: simm,
            ir_dv01_buckets: ir_buckets,
            credit_cs01_buckets: credit_buckets,
            jump_to_default: jtd,
            bermudan_sabr_metric: bermudan_value + sabr_metric,
            sacc_swing_qp_metric: sa_report.ead.val() + swing_value + qp_value,
            slv_rough_worst_of_metric: rough_metric,
            frtb_cvar_execution_metric: execution_metric,
            vix_parisian_robust_metric,
            double_heston_tarf_coco_hrp_metric,
            jy_tba_bergomi_nco_metric,
            cheyette_cmo_qrh_evar_metric,
            asw_mf1f_ah_bess_ep_metric,
            mtmxccs_colva_zabr_weather_enb_metric,
            burnout_cmcds_quintic_afv_cs_schur_metric,
            sinhou_mb_rbergomi_choi_afs_expectile_metric,
            monotone_convex_clark_saccr_rvol_cartea_paleologo_metric: monotone_convex_clark_saccr_rvol_cartea_live_metric(
                &market,
                credit_bump,
            ),
            krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric,
            rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric: self
                .rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_nco_live_metric(),
            mbs_rmv_mo_sig_garch_rie_gp_ierp_metric: self.mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_live_metric(),
            kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric:
                Self::kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_live_metric_from_market(&market, credit_bump),
            mf_rates_wildcard_ctd_pivot_tarf_bergomi_metric: Self::mf_rates_wildcard_ctd_pivot_tarf_bergomi_live_metric_from_market(&market, credit_bump),
            cochrane_ftpl_g2_twist_essvi_offshore_cb_metric: Self::cochrane_ftpl_g2_twist_essvi_offshore_cb_live_metric_from_market(&market, credit_bump),
            veronesi_conv_lando_gen_gatheral_corridor_metric: Self::veronesi_conv_lando_gen_gatheral_corridor_live_metric_from_market_with_bump(
                &market,
                credit_bump,
            ),
            romer_taylor_damodaran_crp_christoffersen_var_metric: Self::romer_taylor_damodaran_crp_christoffersen_var_live_metric_from_market_with_bump(
                &market,
                credit_bump,
            ),
            fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric: Self::fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_live_metric_from_market(&market, credit_bump),
            cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_metric: Self::cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_live_metric_from_market(&market, credit_bump),
            rebonato_cms_leland_toft_ascot_pool_roll_rie_metric: Self::rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric_from_market(&market, credit_bump),
            tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_metric: Self::tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric_from_market(&market, credit_bump),
            hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_metric: Self::hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric_from_market(&market, credit_bump),
            brigo_csa_ctd_callable_oas_jlt_corridor_lsm_metric: Self::brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric_from_market(&market, credit_bump),
            brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_metric: Self::brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric_from_market(&market, credit_bump),
            hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_metric: Self::hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric_from_market(&market, credit_bump),
            brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_metric: Self::brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric_from_market(&market, credit_bump),
            brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_metric: Self::brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric_from_market(&market, credit_bump),
            andersen_cms_special_repo_cds_upfront_ssr_kva_crack_metric: Self::andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric_from_market(&market, credit_bump),
            darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_metric: Self::darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric_from_market(&market, credit_bump),
            henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_metric:
                Self::henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            brigo_diff_swap_callable_muni_bkm_skew_break_cva_metric:
                Self::brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_metric:
                Self::andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_metric:
                Self::litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_metric:
                Self::piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_metric:
                Self::gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_metric:
                Self::kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_metric:
                Self::jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_metric:
                Self::lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_metric:
                Self::jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_metric:
                Self::bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_metric:
                Self::andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_metric:
                Self::piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_metric:
                Self::mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_metric:
                Self::ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_metric:
                Self::smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_metric:
                Self::mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_metric:
                Self::hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_metric:
                Self::mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_metric:
                Self::filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_metric:
                Self::lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_metric:
                Self::cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_metric:
                Self::bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_metric:
                Self::hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_metric:
                Self::klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric_from_market(
                    &market,
                    credit_bump,
                ),
            bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_metric: Self::bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_live_metric_from_market_with_bump(
                &market,
                if self.scenario == ScenarioKind::CreditWideningPlus100Bp {
                    0.01
                } else {
                    0.0
                },
            ),
            blotter_rows,
        }
    }
    fn chart(&self, frame: &WorkstationFrameData) -> BrailleCanvas {
        let mut canvas = BrailleCanvas::new(58, 10);
        let market = self.scenario_market();
        let curves = self.active_curve_system(&market);
        let x32: Vec<f64> = (0..32).map(|i| 0.25 + i as f64 * 9.75 / 31.0).collect();
        match self.chart_mode {
            ChartMode::YieldAndInflationCurves => {
                let (_, inflation) = self.inflation_calibration(self.scenario);
                let ois: Vec<f64> = x32
                    .iter()
                    .map(|&t| curves.curves[0].zero_rate(t).val())
                    .collect();
                let fwd: Vec<f64> = x32
                    .iter()
                    .map(|&t| curves.curves[1].zero_rate(t).val())
                    .collect();
                let cpi: Vec<f64> = x32
                    .iter()
                    .map(|&t| {
                        let month = (((t * 12.0).round() as i32 - 1).rem_euclid(12) + 1) as u8;
                        ZeroCouponInflationSwap {
                            notional: 1.0,
                            fixed_rate: 0.0,
                            maturity_years: t,
                            target_month: Some(month),
                            pay_fixed: true,
                        }
                        .fair_rate(&inflation.curve)
                        .val()
                    })
                    .collect();
                let (lo, hi) = range3(&ois, &fwd, &cpi);
                canvas.plot_series(&x32, &ois, lo, hi, [0, 220, 255]);
                canvas.plot_series(&x32, &fwd, lo, hi, [0, 255, 140]);
                canvas.plot_series(&x32, &cpi, lo, hi, [255, 200, 60]);
            }
            ChartMode::RoughHestonAndFxSmiles => {
                let strikes: Vec<f64> = (0..25).map(|i| 85.0 + i as f64 * 30.0 / 24.0).collect();
                let expiry = 1.0;
                let rate = curves.curves[0].zero_rate(expiry).val();
                let discount = (-rate * expiry).exp();
                let forward = market.eq_spot * (rate * expiry).exp();
                let rough = LiftedHestonParams {
                    v0: market.heston_params.v0,
                    lambda: market.heston_params.kappa,
                    theta: market.heston_params.theta,
                    nu: market.heston_params.sigma_v,
                    rho: market.heston_params.rho,
                    hurst: 0.1,
                    n_factors: 8,
                    r_n: 2.5,
                };
                let heston: Vec<f64> = strikes
                    .iter()
                    .map(|&k| {
                        let price = lifted_heston_cos_price(
                            market.eq_spot,
                            k,
                            expiry,
                            rate,
                            0.0,
                            &rough,
                            64,
                            8.0,
                            16,
                            true,
                        )
                        .val();
                        implied_black_vol(forward, k, expiry, price / discount, true)
                            .map(|v| v.val())
                            .unwrap_or(0.0)
                    })
                    .collect();
                let fx = FxVannaVolgaSmile::new(
                    100.0,
                    0.03,
                    0.01,
                    expiry,
                    FxSmileQuote {
                        atm_vol: 0.1 + 0.1 * market.eq_vol,
                        rr_25d: 0.01,
                        bf_25d: 0.005,
                    },
                );
                let fx_vol: Vec<f64> = strikes.iter().map(|&k| fx.implied_vol(k).val()).collect();
                let (lo, hi) = range2(&heston, &fx_vol);
                canvas.plot_series(&strikes, &heston, lo, hi, [0, 220, 255]);
                canvas.plot_series(&strikes, &fx_vol, lo, hi, [255, 200, 60]);
            }
            ChartMode::XvaExposureProfile => {
                let cube = self.exposure_cube();
                let epe = expected_positive_exposure(&cube);
                let pfe = potential_future_exposure(&cube, 0.975);
                let dim = {
                    let d = self.base_market.ref_date;
                    let active = self.scenario_market();
                    let curve = FlatForwardCurve::new(d, active.ois_quotes[0]);
                    let model = HullWhite1F {
                        a: 0.05,
                        sigma: 0.01 + 0.25 * active.eq_vol,
                    };
                    let swap = NettingSetSwap {
                        notional: 1_000_000.0,
                        fixed_rate: active.ois_quotes[1],
                        payment_times: (1..=5).map(|x| x as f64).collect(),
                        accruals: vec![1.0; 5],
                        is_payer: true,
                    };
                    let times = (1..=16).map(|k| k as f64 * 5.0 / 16.0).collect::<Vec<_>>();
                    ql_xva::compute_chebyshev_dim_and_mva(
                        &model,
                        &curve,
                        &[swap],
                        &times,
                        64,
                        6,
                        0.01,
                        29,
                    )
                    .eim_profile
                };
                let xs: Vec<f64> = cube.times.clone();
                let (lo, hi) = range3(&epe, &pfe, &dim);
                canvas.plot_series(&xs, &epe, lo, hi, [0, 255, 140]);
                canvas.plot_series(&xs, &pfe, lo, hi, [255, 90, 90]);
                canvas.plot_series(&xs, &dim, lo, hi, [180, 130, 255]);
            }
        }
        let _ = frame;
        canvas
    }
    pub fn render_dashboard(&self, frame: &WorkstationFrameData, ansi_color: bool) -> String {
        const LEFT: usize = 57;
        const RIGHT: usize = 58;
        let chart_lines = self.chart(frame).render_lines(ansi_color);
        let mut left_top = Vec::with_capacity(10);
        left_top.push("1. BOOK SUMMARY, XVA & ISDA SIMM".to_string());
        left_top.push(format!(
            "Total NPV: {:>14.2} {}",
            frame.total_npv,
            render_sparkline(&self.npv_history, 16)
        ));
        left_top.push(format!("Intraday PnL: {:>11.2}", frame.intraday_pnl));
        left_top.push(format!(
            "Eq Delta: {:>12.2}  Eq Vega: {:>10.2}",
            frame.eq_delta, frame.eq_vega
        ));
        left_top.push(format!(
            "FX Delta: {:>12.5}  CPI Δ01: {:>10.2}",
            frame.fx_spot_delta, frame.inflation_par_delta_sum
        ));
        left_top.push(format!(
            "XVA CVA/DVA/FVA: {:>7.2}/{:>7.2}/{:>7.2}  MVA:{:>7.2}",
            frame.cva, frame.dva, frame.fva, frame.mva
        ));
        left_top.push(format!("ISDA SIMM IM: {:>14.2}", frame.simm_initial_margin));
        left_top.push(format!("Jump-to-default: {:>12.2}", frame.jump_to_default));
        let cached = frame
            .dag_total_trades
            .saturating_sub(frame.dag_recomputed_trades);
        left_top.push(format!(
            "DAG cache hit: {cached}/{}",
            frame.dag_total_trades
        ));
        left_top.push(format!(
            "Eval: {:>8.0} µs  Mode: {}",
            frame.eval_micros,
            frame.scenario.label()
        ));
        let right_top = chart_lines;

        let max_ir = frame
            .ir_dv01_buckets
            .iter()
            .map(|x| x.1.abs())
            .fold(0.0, f64::max)
            .max(1.0);
        let max_cs = frame
            .credit_cs01_buckets
            .iter()
            .map(|x| x.1.abs())
            .fold(0.0, f64::max)
            .max(1.0);
        let mut left_bottom = Vec::with_capacity(11);
        for (label, value) in &frame.ir_dv01_buckets {
            left_bottom.push(format!(
                "IRS DV01 {label:>3}: {:>9.2} {}",
                value,
                render_diverging_bar(*value, max_ir, 10, ansi_color)
            ));
        }
        for (label, value) in &frame.credit_cs01_buckets {
            left_bottom.push(format!(
                "CDS CS01 {label:>3}: {:>9.2} {}",
                value,
                render_diverging_bar(*value, max_cs, 10, ansi_color)
            ));
        }
        let american = BaroneAdesiWhaleyEngine::price(100.0, 100.0, 0.03, 0.01, 0.2, 1.0, true);
        let erc = equal_risk_contribution_weights(&[vec![0.04, 0.01], vec![0.01, 0.09]]);
        left_bottom.extend([
            format!("19 Bermudan/SABR IFT: {:.2}", frame.bermudan_sabr_metric),
            format!("20 SA-CCR/Swing/QP: {:.2}", frame.sacc_swing_qp_metric),
            format!(
                "21 SLV/rBergomi/WoF: {:.2}",
                frame.slv_rough_worst_of_metric
            ),
            format!("22 FRTB/CVaR/Exec: {:.2}", frame.frtb_cvar_execution_metric),
            format!(
                "23 Monotone/VIX/Parisian/Robust: {:.2}",
                frame.vix_parisian_robust_metric
            ),
            format!(
                "24 Double Heston/TARF/CoCo/HRP: {:.2}",
                frame.double_heston_tarf_coco_hrp_metric
            ),
            format!(
                "25 JY-Inflation/TBA/Bergomi2F/NCO: {:.2}",
                frame.jy_tba_bergomi_nco_metric
            ),
            format!(
                "26 Cheyette/CMO/QRH/BaseCorr/EVaR: {:.2}",
                frame.cheyette_cmo_qrh_evar_metric
            ),
            format!(
                "27 ASW/MF1F/AH-LV/BESS/Idzorek-DP: {:.2}",
                frame.asw_mf1f_ah_bess_ep_metric
            ),
            format!(
                "28 MtM-XCCS/COLVA/ZABR/Weather/ENB: {:.2}",
                frame.mtmxccs_colva_zabr_weather_enb_metric
            ),
            format!(
                "29 Burnout/CMCDS/Quintic/AFV/Schur: {:.2}",
                frame.burnout_cmcds_quintic_afv_cs_schur_metric
            ),
            format!(
                "30 SinhOU/MB/mrBergomi/Choi/EXVaR: {:.2}",
                frame.sinhou_mb_rbergomi_choi_afs_expectile_metric
            ),
            format!(
                "31 HW-MC/Clark/SACCR/rVol/Cartea/Paleo: {:.2}",
                frame.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
            ),
            format!(
                "32 KRD/LPI/CoCo/Bergomi2F/CGMY/Guéant/Meucci: {:.2}",
                frame.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
            ),
            format!(
                "33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML: {:.2}",
                frame.rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric
            ),
            format!(
                "34 MBS/RMV/MO/CB/Sig/HN/POT/RIE/PIN/GP/BC/IERP: {:.2}",
                frame.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric
            ),
            format!(
                "35 Kalman/VECM/HJ/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin: {:.2}",
                frame
                    .kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric
            ),
            format!("36 MF-Rates/Wildcard/TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC: {:.2}", frame.mf_rates_wildcard_ctd_pivot_tarf_bergomi_metric),
            format!("37 FTPL/G2/eSSVI/CB/CoCo/MPOR/Double-Heston/Curran/Warin/Roll/Busseti/KKT: {:.2}", frame.cochrane_ftpl_g2_twist_essvi_offshore_cb_metric),
            format!("38 IRRBB/RFR/XCCY/BG/CoCo/SIMM/rBergomi/Spread/BESS/Dark/Risk/Diff-PCA: {:.2}", frame.bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_metric),
            format!("39 Veronesi/Lando/Corridor/Russian/Bates/Koziol/MFG/BCBS/CJP/ETF/IC/IS: {:.2}", frame.veronesi_conv_lando_gen_gatheral_corridor_metric),
            format!("40 Taylor/Damodaran/VaR/LVaR/CGO/Regret/Beneish/Campisi/Heckman/Elastic/ResetCB/FSB: {:.2}", frame.romer_taylor_damodaran_crp_christoffersen_var_metric),
            format!("41 KRD-LSC/CMO-PAC/Asquith-CB/ASU2024/ASC470/ColVA/Hill/Bouchaud/Karnosky/DGTW/Brock/Penman: {:.2}", frame.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric),
            format!("42 Cheyette/NDIRS/ZABR/QRH/CMCDS/WWR/HHW/CB/EVaR/Schur/BEKK/Quality: {:.2}", frame.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_metric),
            format!("43 Rebonato-CMS/Leland-Toft/Calamos/RIE/Arrow-Debreu: {:.2}", frame.rebonato_cms_leland_toft_ascot_pool_roll_rie_metric),
            format!("44 Tuckman-Repo/BK-FVA/DECS/MBS-IO-PO/Kyle-PIN/GP-Aim/GJR-DCC/LPPL/Spectral-Risk/Brinson-Carino/Prelec-CPT/Path-Sig: {:.2}", frame.tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_metric),
            format!("45 Hagan-West-MC/CTD-Switch/AFV-CB/Kth-Default/SA-CVA/GS-Kirk/FX-VV/Napoleon/TC-FLAM/Meucci/KRR/EVA: {:.2}", frame.hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_metric),
            format!("46 Multi-CSA-CTD/Callable-OAS/JLT-Migration/Corridor-VarSwap/LSM-Dual/Cartea-LOB-MM/Lo-MacKinlay-VR/Lopez-DSR/Roncalli-ERC/Ilmanen-VMC/Penman-ReOI/McCall: {:.2}", frame.brigo_csa_ctd_callable_oas_jlt_corridor_lsm_metric),
            format!("47 NSS-KRD/TIPS-Floor/RMV-RT-RFV/Sticky-LV/SIMM-CCP/Schwartz-Smith/FX-DNT/Timer/BL-Woodbury/AC-Sinh/ENet/NK: {:.2}", frame.brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_metric),
            format!("48 Tension-Curve/Duffie-Lando/PAC-CMO/FM-Shanken/aSVI/BGK-Barrier/Cross-Impact/Karnosky/Merton/CV/Ohlson/Revolve: {:.2}", frame.hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_metric),
            format!("49 XCCY-MtM/Spec-Pool/Marshall-Olkin/VIX-VVIX/WWR-CVA/Gibson-Schwartz/FX-FVA/Chooser/FLAM/ES-Euler/SURE/Habit: {:.2}", frame.brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_metric),
            format!("50 OIS-Fallback/FRN/CDS-Swaption/Bergomi-2F/SA-CCR/Spark-Spread/PRDC/Geske/Alpha-Decay/ENB/VECM/Bansal-Yaron: {:.2}", frame.brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_metric),
            format!("51 CMS-Rep/Spec-Repo/CDS-Upfront/SSR-Tanaka/SA-CCR-KVA/Crack/PRDC-Cap/ESO/Spinu-ERC/CVaR/PLS-Krylov/FOMC: {:.2}", frame.andersen_cms_special_repo_cds_upfront_ssr_kva_crack_metric),
            format!("52 Step-OIS-Turn/CB-SoftCall/CDX-Skew/VarSwap-Trunc/FVA-Overlap/Working-Storage/FX-Fly/Lookback/Ledoit-Wolf/Almgren-CARA/ESTAR/ACM-CP: {:.2}", frame.darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_metric),
            format!("53 Tenor-Basis/MSR-Recap/CDO-BaseCorr/Fwd-Smile/CSA-MTA/Gas-Swing/FX-KIKO/Stulz-Rainbow/Jagannathan-Ma/Obizhaeva-Wang/Newey-West/Gabaix-Disaster: {:.2}", frame.henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_metric),
            format!("54 Diff-Swap/Muni-YTW/Step-Up-Bond/BKM-SKEW/Break-CVA/Power-Skew/FX-Fader/Extendible/Kan-Zhou/Sharpe-RBSA/DM-HLN/Laubach-r*: {:.2}", frame.brigo_diff_swap_callable_muni_bkm_skew_break_cva_metric),
            format!("55 SOFR-Fut/Hayre-MBS/CDX-Opt/Dupire-LV/SIMM-MVA/RSS-Stockout/CIP-Basis/Carr-Chou/Pastor-Stambaugh/Treynor-Black/Engle-Granger/Epstein-Zin: {:.2}", frame.andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_metric),
            format!("56 PCA-Fly/TF-CB/Double-t-CDO/Lee-Wings/CDS-Contagion/Schwartz-OU/FX-DCD/Passport/GRS-FTest/Sqrt-Impact/BNS-HAR-CJ/LL-cay: {:.2}", frame.litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_metric),
            format!("57 Piterbarg-Opt/Tuckman-Roll/Merton-KMV/Bartlett-SABR/Pykhtin-MPoR/Deng-FTR/FX-Digital/Conze-LB/Hot-Spots-ES/Frazzini-BAB/PT-Test/Lucas-MP: {:.2}", frame.piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_metric),
            format!("58 GSW-Svensson/MBS-ESF/Vasicek-IRB-GA/Carr-Wu-Vol/Marginal-CVA/Working-Storage/FX-Fwd-Extra/Exch-on-Exch/HLZ-FDR/HS-Kyle/GARCH-Kurt/CS-CAPE: {:.2}", frame.gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_metric),
            format!("59 Kim-Wright-Shadow/Brennan-Schwartz-CB/Duffie-RMV/Fukasawa-Gamma/Burgard-Kjaer-FVA/CD-Spark/FX-TARF/BGK-Lookback/Novy-Marx-FF5/Acerbi-SRM/Corsi-HAR-RV/BLL-News-Noise: {:.2}", frame.kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_metric),
            format!("60 JPS-Unspanned-Macro/Kalotay-Refunding/Hull-White-NtD/Derman-Sticky/HW-WWR-CVA/GS-Samuelson/Reiner-Quanto-FX/TW-Asian/AMP-VME/GP-Aim/Bai-Ng-ICp2/BGG-Accelerator: {:.2}", frame.jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_metric),
            format!("61 Lyashenko-FMM-Caplet/Tuckman-CPI-Seasonality/Brigo-Alfonsi-CIR/Fouque-SV/Gregory-FtD/Eydeland-VPP/Wystup-Corridor/RR-Barrier/LNS-GLS-R2/Bellini-Expectile/BNS-Quarticity/Kiyotaki-Moore: {:.2}", frame.lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_metric),
            format!("62 JSZ-Canonical-GDTSM/Tsiveriotis-Fernandes-CB/Laurent-Gregory-CDO/Broadie-Jain-VarSwap/Albanese-SA-CCR-KVA/RSS-Storage/DTV-CIP/Wilmott-Cliquet/BSV-PPP/Ahmadi-Javid-EVaR/Engle-Kroner-BEKK/Gabaix-Disaster: {:.2}", frame.jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_metric),
            format!("63 Bianchetti-Two-Curve-Quanto/BDT-Callable/Schonbucher-Swaption/Bates-SVJ/Brigo-Pallavicini-CCP/DMT-Crack/Piterbarg-PRDC/BGK-Barrier/Connor-Korajczyk-APT/AFS-Impact/Hamilton-MS/Bansal-Yaron: {:.2}", frame.bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_metric),
            format!("64 Andersen-Piterbarg-CMS-Spread/Fabozzi-Z-Bond/Duffie-Garleanu-CDO/Carr-Lee-VolSwap/Green-LR-KVA/Cortazar-Schwartz-3F/Fujii-XCCS/Pechtl-Haug-TimeSwitch/Boudt-CF-MVaR/Rockafellar-Quadrangle/Hasbrouck-IS/He-Krishnamurthy: {:.2}", frame.andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_metric),
            format!("65 Piterbarg-CTD-Collateral/Tuckman-Make-Whole/Brigo-Morini-CMCDS/Hagan-Shifted-SABR/Green-Kenyon-TVA/Eydeland-Interconnector/Wystup-Clark-Strangle/Kunitomo-Ikeda-Barrier/CLMX-Var/Kyle-Lambda/Hansen-Lunde-RV/Goldstein-Pauzner: {:.2}", frame.piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_metric),
            format!("66 Mercurio-FVA-Swap/Ingersoll-CB/Hull-White-Implied-Copula/Gatheral-Jacquier-SVI-JW/Burgard-Kjaer-XVA/Benth-HDD-CDD/Baba-Packer-FX-Basis/Conze-Viswanathan-Lookback/He-Litterman-BL/GHS-IS-CVaR/Engle-Granger-Johansen/Gabaix-Disasters: {:.2}", frame.mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_metric),
            format!("67 Ametrano-Bianchetti-Tenor-Basis/Fabozzi-Hayre-IO-PO/Hull-Brigo-WWR-CDS/Bossens-Vanna-Volga-Barrier/Gregory-SIMM-MVA/Carmona-Clean-Dark-Spread/Reiner-Haug-Quanto/Hart-Ross-Soft-Barrier/Pastor-Stambaugh-Alpha/Obizhaeva-Wang-LOB/DMW-HLN/Piazzesi-Schneider: {:.2}", frame.ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_metric),
            format!("68 Smith-Wilson-UFR/Kalotay-Williams-Fabozzi/Schonbucher-Clayton-CDO/Buehler-Variance-Curve/Pykhtin-Rosen-Marginal-CVA/Hilliard-Reis-Futures/Gabaix-Maggiori-FX/Cheuk-Vorst-Shout/Grossman-Zhou-Kelly/GLS-Tasche-RAROC/Andrews-Ploberger/Jermann-Quadrini: {:.2}", frame.smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_metric),
            format!("69 McCulloch-Tax-Curve/Andersen-Buffum-CB/OKane-Recovery-Lock/Derman-Kani-Tree/Kenyon-Green-LVA/Boogert-de-Jong-Gas/Garman-Kohlhagen-BAW/Hui-Haug-Binary/Kan-Robotti-HJ/Cartea-Jaimungal-OF/Ait-Sahalia-Jacod/Wachter-Disaster: {:.2}", frame.mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_metric),
            format!("70 Hagan-West-Monotone-Convex/Fabozzi-Effective-Duration-OAS/Frye-Altman-LGD/Gatheral-Jacquier-SVI/Gregory-Brigo-BCVA/Brennan-Schwartz-Mine/Wystup-FX-Fader/Broadie-Detemple-Cap/Meucci-Entropy/Cont-Deguest-Scandolo/Pesaran-Shin-Smith/CEE-NK: {:.2}", frame.hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_metric),
            format!("71 Mercurio-Russo-IBOR-OIS-Basis/Jarrow-Yildirim-YoY-Inflation/Duffie-Lando-CDS/Bayer-Friz-Gatheral-Skew/Andersen-Pykhtin-Sokol-MPoR/Bjerksund-Stensland-Crack/Carr-Chou-PCS/Chesney-Jeanblanc-Yor-Parisian/DeMiguel-Brodie-Norm/Gueant-Lehalle-Tapia-MM/Hansen-White-SPA/Brunnermeier-Sannikov-Paradox: {:.2}", frame.mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_metric),
            format!("72 Filipovic-Larsson-Trolle-LRTSM/Duffie-Krishnamurthy-Repo/Duffie-Saita-Wang-Default/Carr-Madan-Neuberger-LogContract/Hull-White-Burgard-Kjaer-FVA/Casassus-Collin-Dufresne-CY/Lustig-Roussanov-Verdelhan-FX/Carr-Lee-Corridor-VarSwap/Garlappi-Uppal-Wang-MultiPrior/McNeil-Frey-EVT-GARCH/Patton-SJC-Copula/Pastor-Veronesi-Political: {:.2}", frame.filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_metric),
            format!("73 Lyashenko-Mercurio-RFR/Pennacchi-Sundaresan-CoCo/Pedersen-Brigo-CDX-FEP/Broadie-Jain-Sepp-SVJ-VolSwap/Brigo-Capponi-Pallavicini-WWR/Routledge-Seppi-Spatt-Storage/Piterbarg-PRDC-Hybrid/Broadie-Glasserman-Kou-Barrier/Brandt-Santa-Clara-Valkanov-PPP/Gatheral-Schied-Slynko-Impact/Corsi-BPQ-HARQ/Lettau-Ludvigson-CAY: {:.2}", frame.lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_metric),
            format!("74 Cieslak-Povala-Cycle/Gabaix-Krishnamurthy-Vigneron-MBS/Errais-Giesecke-Goldberg-Hawkes/Martin-Bondarenko-SVIX/Crepey-Bichuch-Capponi-BSDE/Trolle-Schwartz-USV/Ivashina-Scharfstein-Stein-CIP/Wilmott-Windcliff-Cliquet/Ledoit-Wolf-Nonlinear-Shrinkage/Adrian-Brunnermeier-CoVaR/Bollerslev-Tauchen-Zhou-VRP/Campbell-Vuolteenaho-Beta: {:.2}", frame.cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_metric),
            format!("75 Bianchetti-Mercurio-FRA-Quanto/Acharya-Amihud-Bharath-Liq/Brigo-Alfonsi-SSRD-CDS/Fukasawa-Skew-VarSwap/Andersen-Duffie-Song-FVA/Pirrong-Inventory-Spread/Lipton-Castagna-FX-DNT/Linetsky-Step-Barrier/Koijen-Yogo-Demand-System/Kyle-Obizhaeva-Invariance/Gu-Kelly-Xiu-OOS-R2/Drechsler-Savov-Schnabl-Deposit: {:.2}", frame.bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_metric),
            format!("76 Hanson-Stein-Duration/Fleckenstein-Longstaff-Lustig-TIPS/Collin-Dufresne-Goldstein-Yang-CDO/Todorov-Spot-Vol/Albanese-Caenazzo-Crepey-KVA/Acharya-Lochstoer-Ramadorai-Hedging/Brunnermeier-Nagel-Pedersen-FX/Broadie-Yamamoto-FGT-Asian/Kelly-Pruitt-Su-IPCA/Brownlees-Engle-SRISK/ABDL-ARFIMA-RV/Farhi-Gabaix-Disaster-FX: {:.2}", frame.hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_metric),
            format!("77 Klingler-Syrstad-SOFR-Spike/Longstaff-Mithal-Neis-CDS-Bond/Pan-Singleton-Sovereign-CDS/Medvedev-Scaillet-Short-Smile/Duffie-Zhu-CCP-Netting/Gorton-Hayashi-Rouwenhorst-Inventory/Della-Corte-Ramadorai-Sarno-FX-VRP/Pelsser-Laplace-Double-Barrier/Moreira-Muir-Vol-Managed/Budish-Cramton-Shim-FBA/AMZ-Two-Scales-RV/Gertler-Karadi-QE: {:.2}", frame.klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_metric),
            format!(
                "DAG recomputed: {}/{}",
                frame.dag_recomputed_trades, frame.dag_total_trades
            ),
            format!("Chart: {}", self.chart_mode.label()),
            format!("American / Hybrid: {:.4}", american),
            format!("BL / Risk Parity ERC: {:.3}/{:.3}", erc[0], erc[1]),
        ]);
        // There are six IR buckets, five credit buckets and ten telemetry
        // lines (including Tab 23).  Render all 21 rows instead of silently
        // dropping the final BL/Risk-Parity line.

        let mut right_bottom = vec![];
        right_bottom.extend(
            frame
                .blotter_rows
                .iter()
                .map(|(name, value, tag)| format!("{name:<20} {:>12.2} {tag:<10}", value)),
        );
        right_bottom.push("32 KRD/LPI/CoCo/Bergomi/CGMY/AAD live".to_string());
        right_bottom.push("33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML live".to_string());
        right_bottom.push("34 Davidson-Levin MBS IO/PO/Duffie-Singleton RMV/MO-XVA/Sirbu-Shreve CB/PathSig/HN-DCC-GARCH/POT-GPD/RIE/PIN/GP-Aim/Brinson/IERP".to_string());
        right_bottom.push(format!(
            "35 Kalman/VECM/HJ/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin: {:.4}",
            self.kalman_rts_vecm_hj_sdf_ifrs9_cpt_hasbrouck_live_metric()
        ));
        right_bottom.push(format!("36 MF-Rates/Wildcard/TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC: {:.4}", self.mf_rates_wildcard_ctd_pivot_tarf_bergomi_live_metric()));
        right_bottom.push(format!(
            "37 FTPL/G2/eSSVI/CB/CoCo/MPOR/Double-Heston/Curran/Warin/Roll/Busseti/KKT: {:.4}",
            self.cochrane_ftpl_g2_twist_essvi_offshore_cb_live_metric()
        ));
        right_bottom.push(format!(
            "38 IRRBB/RFR/XCCY/BG/CoCo/SIMM/rBergomi/Spread/BESS/Dark/Risk/Diff-PCA: {:.4}",
            self.bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_live_metric()
        ));
        right_bottom.push(format!(
            "39 Veronesi/Lando/Corridor/Russian/Bates/Koziol/MFG/BCBS/CJP/ETF/IC/IS: {:.4}",
            self.veronesi_conv_lando_gen_gatheral_corridor_live_metric()
        ));
        right_bottom.push(format!("40 Taylor/Damodaran/VaR/LVaR/CGO/Regret/Beneish/Campisi/Heckman/Elastic/ResetCB/FSB: {:.4}", self.romer_taylor_damodaran_crp_christoffersen_var_live_metric()));
        right_bottom.push(format!("41 KRD-LSC/CMO-PAC/Asquith-CB/ASU2024/ASC470/ColVA/Hill/Bouchaud/Karnosky/DGTW/Brock/Penman: {:.4}", self.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_live_metric()));
        right_bottom.push(format!(
            "42 Cheyette/NDIRS/ZABR/QRH/CMCDS/WWR/HHW/CB/EVaR/Schur/BEKK/Quality: {:.4}",
            self.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_live_metric()
        ));
        right_bottom.push(format!(
            "43 Rebonato-CMS/Leland-Toft/Calamos/RIE/Arrow-Debreu: {:.4}",
            self.rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric()
        ));
        right_bottom.push(format!("44 Tuckman-Repo/BK-FVA/DECS/MBS-IO-PO/Kyle-PIN/GP-Aim/GJR-DCC/LPPL/Spectral-Risk/Brinson-Carino/Prelec-CPT/Path-Sig: {:.4}", self.tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric()));
        right_bottom.push(format!("45 Hagan-West-MC/CTD-Switch/AFV-CB/Kth-Default/SA-CVA/GS-Kirk/FX-VV/Napoleon/TC-FLAM/Meucci/KRR/EVA: {:.4}", self.hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric()));
        right_bottom.push(format!("46 Multi-CSA-CTD/Callable-OAS/JLT-Migration/Corridor-VarSwap/LSM-Dual/Cartea-LOB-MM/Lo-MacKinlay-VR/Lopez-DSR/Roncalli-ERC/Ilmanen-VMC/Penman-ReOI/McCall: {:.4}", self.brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric()));
        right_bottom.push(format!("47 NSS-KRD/TIPS-Floor/RMV-RT-RFV/Sticky-LV/SIMM-CCP/Schwartz-Smith/FX-DNT/Timer/BL-Woodbury/AC-Sinh/ENet/NK: {:.4}", self.brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric()));
        right_bottom.push(format!("48 Tension-Curve/Duffie-Lando/PAC-CMO/FM-Shanken/aSVI/BGK-Barrier/Cross-Impact/Karnosky/Merton/CV/Ohlson/Revolve: {:.4}", self.hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric()));
        right_bottom.push(format!("49 XCCY-MtM/Spec-Pool/Marshall-Olkin/VIX-VVIX/WWR-CVA/Gibson-Schwartz/FX-FVA/Chooser/FLAM/ES-Euler/SURE/Habit: {:.4}", self.brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric()));
        right_bottom.push(format!("50 OIS-Fallback/FRN/CDS-Swaption/Bergomi-2F/SA-CCR/Spark-Spread/PRDC/Geske/Alpha-Decay/ENB/VECM/Bansal-Yaron: {:.4}", self.brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric()));
        right_bottom.push(format!("51 CMS-Rep/Spec-Repo/CDS-Upfront/SSR-Tanaka/SA-CCR-KVA/Crack/PRDC-Cap/ESO/Spinu-ERC/CVaR/PLS-Krylov/FOMC: {:.4}", self.andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric()));
        right_bottom.push(format!("52 Step-OIS-Turn/CB-SoftCall/CDX-Skew/VarSwap-Trunc/FVA-Overlap/Working-Storage/FX-Fly/Lookback/Ledoit-Wolf/Almgren-CARA/ESTAR/ACM-CP: {:.4}", self.darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric()));
        right_bottom.push(format!("53 Tenor-Basis/MSR-Recap/CDO-BaseCorr/Fwd-Smile/CSA-MTA/Gas-Swing/FX-KIKO/Stulz-Rainbow/Jagannathan-Ma/Obizhaeva-Wang/Newey-West/Gabaix-Disaster: {:.4}", self.henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric()));
        right_bottom.push(format!("54 Diff-Swap/Muni-YTW/Step-Up-Bond/BKM-SKEW/Break-CVA/Power-Skew/FX-Fader/Extendible/Kan-Zhou/Sharpe-RBSA/DM-HLN/Laubach-r*: {:.4}", self.brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric()));
        right_bottom.push(format!("55 SOFR-Fut/Hayre-MBS/CDX-Opt/Dupire-LV/SIMM-MVA/RSS-Stockout/CIP-Basis/Carr-Chou/Pastor-Stambaugh/Treynor-Black/Engle-Granger/Epstein-Zin: {:.4}", self.andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric()));
        right_bottom.push(format!("56 PCA-Fly/TF-CB/Double-t-CDO/Lee-Wings/CDS-Contagion/Schwartz-OU/FX-DCD/Passport/GRS-FTest/Sqrt-Impact/BNS-HAR-CJ/LL-cay: {:.4}", self.litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric()));
        right_bottom.push(format!("57 Piterbarg-Opt/Tuckman-Roll/Merton-KMV/Bartlett-SABR/Pykhtin-MPoR/Deng-FTR/FX-Digital/Conze-LB/Hot-Spots-ES/Frazzini-BAB/PT-Test/Lucas-MP: {:.4}", self.piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric()));
        right_bottom.push(format!("58 GSW-Svensson/MBS-ESF/Vasicek-IRB-GA/Carr-Wu-Vol/Marginal-CVA/Working-Storage/FX-Fwd-Extra/Exch-on-Exch/HLZ-FDR/HS-Kyle/GARCH-Kurt/CS-CAPE: {:.4}", self.gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric()));
        right_bottom.push(format!("59 Kim-Wright-Shadow/Brennan-Schwartz-CB/Duffie-RMV/Fukasawa-Gamma/Burgard-Kjaer-FVA/CD-Spark/FX-TARF/BGK-Lookback/Novy-Marx-FF5/Acerbi-SRM/Corsi-HAR-RV/BLL-News-Noise: {:.4}", self.kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric()));
        right_bottom.push(format!("60 JPS-Unspanned-Macro/Kalotay-Refunding/Hull-White-NtD/Derman-Sticky/HW-WWR-CVA/GS-Samuelson/Reiner-Quanto-FX/TW-Asian/AMP-VME/GP-Aim/Bai-Ng-ICp2/BGG-Accelerator: {:.4}", self.jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric()));
        right_bottom.push(format!("61 Lyashenko-FMM-Caplet/Tuckman-CPI-Seasonality/Brigo-Alfonsi-CIR/Fouque-SV/Gregory-FtD/Eydeland-VPP/Wystup-Corridor/RR-Barrier/LNS-GLS-R2/Bellini-Expectile/BNS-Quarticity/Kiyotaki-Moore: {:.4}", self.lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric()));
        right_bottom.push(format!("62 JSZ-Canonical-GDTSM/Tsiveriotis-Fernandes-CB/Laurent-Gregory-CDO/Broadie-Jain-VarSwap/Albanese-SA-CCR-KVA/RSS-Storage/DTV-CIP/Wilmott-Cliquet/BSV-PPP/Ahmadi-Javid-EVaR/Engle-Kroner-BEKK/Gabaix-Disaster: {:.4}", self.jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric()));
        right_bottom.push(format!("63 Bianchetti-Two-Curve-Quanto/BDT-Callable/Schonbucher-Swaption/Bates-SVJ/Brigo-Pallavicini-CCP/DMT-Crack/Piterbarg-PRDC/BGK-Barrier/Connor-Korajczyk-APT/AFS-Impact/Hamilton-MS/Bansal-Yaron: {:.4}", self.bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric()));
        right_bottom.push(format!("64 Andersen-Piterbarg-CMS-Spread/Fabozzi-Z-Bond/Duffie-Garleanu-CDO/Carr-Lee-VolSwap/Green-LR-KVA/Cortazar-Schwartz-3F/Fujii-XCCS/Pechtl-Haug-TimeSwitch/Boudt-CF-MVaR/Rockafellar-Quadrangle/Hasbrouck-IS/He-Krishnamurthy: {:.4}", self.andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric()));
        right_bottom.push(format!("65 Piterbarg-CTD-Collateral/Tuckman-Make-Whole/Brigo-Morini-CMCDS/Hagan-Shifted-SABR/Green-Kenyon-TVA/Eydeland-Interconnector/Wystup-Clark-Strangle/Kunitomo-Ikeda-Barrier/CLMX-Var/Kyle-Lambda/Hansen-Lunde-RV/Goldstein-Pauzner: {:.4}", self.piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric()));
        right_bottom.push(format!("66 Mercurio-FVA-Swap/Ingersoll-CB/Hull-White-Implied-Copula/Gatheral-Jacquier-SVI-JW/Burgard-Kjaer-XVA/Benth-HDD-CDD/Baba-Packer-FX-Basis/Conze-Viswanathan-Lookback/He-Litterman-BL/GHS-IS-CVaR/Engle-Granger-Johansen/Gabaix-Disasters: {:.4}", self.mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric()));
        right_bottom.push(format!("67 Ametrano-Bianchetti-Tenor-Basis/Fabozzi-Hayre-IO-PO/Hull-Brigo-WWR-CDS/Bossens-Vanna-Volga-Barrier/Gregory-SIMM-MVA/Carmona-Clean-Dark-Spread/Reiner-Haug-Quanto/Hart-Ross-Soft-Barrier/Pastor-Stambaugh-Alpha/Obizhaeva-Wang-LOB/DMW-HLN/Piazzesi-Schneider: {:.4}", self.ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric()));
        right_bottom.push(format!("68 Smith-Wilson-UFR/Kalotay-Williams-Fabozzi/Schonbucher-Clayton-CDO/Buehler-Variance-Curve/Pykhtin-Rosen-Marginal-CVA/Hilliard-Reis-Futures/Gabaix-Maggiori-FX/Cheuk-Vorst-Shout/Grossman-Zhou-Kelly/GLS-Tasche-RAROC/Andrews-Ploberger/Jermann-Quadrini: {:.4}", self.smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric()));
        right_bottom.push(format!("69 McCulloch-Tax-Curve/Andersen-Buffum-CB/OKane-Recovery-Lock/Derman-Kani-Tree/Kenyon-Green-LVA/Boogert-de-Jong-Gas/Garman-Kohlhagen-BAW/Hui-Haug-Binary/Kan-Robotti-HJ/Cartea-Jaimungal-OF/Ait-Sahalia-Jacod/Wachter-Disaster: {:.4}", self.mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric()));
        right_bottom.push(format!("70 Hagan-West-Monotone-Convex/Fabozzi-Effective-Duration-OAS/Frye-Altman-LGD/Gatheral-Jacquier-SVI/Gregory-Brigo-BCVA/Brennan-Schwartz-Mine/Wystup-FX-Fader/Broadie-Detemple-Cap/Meucci-Entropy/Cont-Deguest-Scandolo/Pesaran-Shin-Smith/CEE-NK: {:.4}", self.hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric()));
        right_bottom.push(format!("71 Mercurio-Russo-IBOR-OIS-Basis/Jarrow-Yildirim-YoY-Inflation/Duffie-Lando-CDS/Bayer-Friz-Gatheral-Skew/Andersen-Pykhtin-Sokol-MPoR/Bjerksund-Stensland-Crack/Carr-Chou-PCS/Chesney-Jeanblanc-Yor-Parisian/DeMiguel-Brodie-Norm/Gueant-Lehalle-Tapia-MM/Hansen-White-SPA/Brunnermeier-Sannikov-Paradox: {:.4}", self.mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric()));
        right_bottom.push(format!("72 Filipovic-Larsson-Trolle-LRTSM/Duffie-Krishnamurthy-Repo/Duffie-Saita-Wang-Default/Carr-Madan-Neuberger-LogContract/Hull-White-Burgard-Kjaer-FVA/Casassus-Collin-Dufresne-CY/Lustig-Roussanov-Verdelhan-FX/Carr-Lee-Corridor-VarSwap/Garlappi-Uppal-Wang-MultiPrior/McNeil-Frey-EVT-GARCH/Patton-SJC-Copula/Pastor-Veronesi-Political: {:.4}", self.filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric()));
        right_bottom.push(format!("73 Lyashenko-Mercurio-RFR/Pennacchi-Sundaresan-CoCo/Pedersen-Brigo-CDX-FEP/Broadie-Jain-Sepp-SVJ-VolSwap/Brigo-Capponi-Pallavicini-WWR/Routledge-Seppi-Spatt-Storage/Piterbarg-PRDC-Hybrid/Broadie-Glasserman-Kou-Barrier/Brandt-Santa-Clara-Valkanov-PPP/Gatheral-Schied-Slynko-Impact/Corsi-BPQ-HARQ/Lettau-Ludvigson-CAY: {:.4}", self.lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric()));
        right_bottom.push(format!("74 Cieslak-Povala-Cycle/Gabaix-Krishnamurthy-Vigneron-MBS/Errais-Giesecke-Goldberg-Hawkes/Martin-Bondarenko-SVIX/Crepey-Bichuch-Capponi-BSDE/Trolle-Schwartz-USV/Ivashina-Scharfstein-Stein-CIP/Wilmott-Windcliff-Cliquet/Ledoit-Wolf-Nonlinear-Shrinkage/Adrian-Brunnermeier-CoVaR/Bollerslev-Tauchen-Zhou-VRP/Campbell-Vuolteenaho-Beta: {:.4}", self.cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric()));
        right_bottom.push(format!("75 Bianchetti-Mercurio-FRA-Quanto/Acharya-Amihud-Bharath-Liq/Brigo-Alfonsi-SSRD-CDS/Fukasawa-Skew-VarSwap/Andersen-Duffie-Song-FVA/Pirrong-Inventory-Spread/Lipton-Castagna-FX-DNT/Linetsky-Step-Barrier/Koijen-Yogo-Demand-System/Kyle-Obizhaeva-Invariance/Gu-Kelly-Xiu-OOS-R2/Drechsler-Savov-Schnabl-Deposit: {:.4}", self.bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric()));
        right_bottom.push(format!("76 Hanson-Stein-Duration/Fleckenstein-Longstaff-Lustig-TIPS/Collin-Dufresne-Goldstein-Yang-CDO/Todorov-Spot-Vol/Albanese-Caenazzo-Crepey-KVA/Acharya-Lochstoer-Ramadorai-Hedging/Brunnermeier-Nagel-Pedersen-FX/Broadie-Yamamoto-FGT-Asian/Kelly-Pruitt-Su-IPCA/Brownlees-Engle-SRISK/ABDL-ARFIMA-RV/Farhi-Gabaix-Disaster-FX: {:.4}", self.hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric()));
        right_bottom.push(format!("77 Klingler-Syrstad-SOFR-Spike/Longstaff-Mithal-Neis-CDS-Bond/Pan-Singleton-Sovereign-CDS/Medvedev-Scaillet-Short-Smile/Duffie-Zhu-CCP-Netting/Gorton-Hayashi-Rouwenhorst-Inventory/Della-Corte-Ramadorai-Sarno-FX-VRP/Pelsser-Laplace-Double-Barrier/Moreira-Muir-Vol-Managed/Budish-Cramton-Shim-FBA/AMZ-Two-Scales-RV/Gertler-Karadi-QE: {:.4}", self.klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric()));
        let bottom_rows = left_bottom.len().max(right_bottom.len());
        left_bottom.resize(bottom_rows, String::new());
        right_bottom.resize(bottom_rows, String::new());

        let mut out = String::new();
        let top_left = pad_visible(
            &format!("─ ql-rs QUANT RISK WORKSTATION  Tick #{:04} ", frame.tick),
            LEFT,
        );
        let top_right = pad_visible(
            &format!("─ 2. BRAILLE CANVAS: {} ", self.chart_mode.label()),
            RIGHT,
        );
        let _ = writeln!(out, "╭{}┬{}╮", top_left, top_right);
        for i in 0..10 {
            let l = left_top.get(i).map(String::as_str).unwrap_or("");
            let r = right_top.get(i).map(String::as_str).unwrap_or("");
            let _ = writeln!(out, "│{}│{}│", pad_visible(l, LEFT), pad_visible(r, RIGHT));
        }
        let _ = writeln!(
            out,
            "├{}┼{}┤",
            pad_visible(" 3. IFT BUCKETED DV01 & CREDIT CS01 LADDERS", LEFT),
            pad_visible(" 4. CROSS-ASSET LIVE BLOTTER & SOTA TELEMETRY", RIGHT)
        );
        for i in 0..bottom_rows {
            let _ = writeln!(
                out,
                "│{}│{}│",
                pad_visible(&left_bottom[i], LEFT),
                pad_visible(&right_bottom[i], RIGHT)
            );
        }
        let _ = writeln!(out, "╰{}┴{}╯", "─".repeat(LEFT), "─".repeat(RIGHT));
        out
    }
}
impl Default for DeskWorkstation {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod domain_suite_26_tests {
    use super::*;
    #[test]
    fn suite26_tab_is_registered() {
        assert_eq!(TabId::ALL.len(), 77);
        assert!(TabId::ALL.contains(&TabId::DoubleHestonTarfCocoHrp));
        assert!(TabId::ALL.contains(&TabId::JyInflationTbaBergomiNco));
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert!(frame.double_heston_tarf_coco_hrp_metric > 100.0);
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("24 Double Heston/TARF/CoCo/HRP")
        );
    }
}

#[cfg(test)]
mod domain_suite_27_tests {
    use super::*;
    #[test]
    fn suite27_tab_and_metric_are_live() {
        assert!(TabId::ALL.contains(&TabId::JyInflationTbaBergomiNco));
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert!(frame.jy_tba_bergomi_nco_metric.is_finite());
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("25 JY-Inflation/TBA/Bergomi2F/NCO")
        );
        ws.set_scenario(ScenarioKind::RatesPlus50Bp);
        let stressed = ws.step_tick();
        assert_ne!(
            stressed.jy_tba_bergomi_nco_metric,
            frame.jy_tba_bergomi_nco_metric
        );
    }
}

#[cfg(test)]
mod domain_suite_28_tests {
    use super::*;

    #[test]
    fn suite28_tab_and_metric_are_live() {
        assert_eq!(TabId::ALL.len(), 77);
        assert!(TabId::ALL.contains(&TabId::CheyetteCmoQrhBaseCorrCdar));
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert!(frame.cheyette_cmo_qrh_evar_metric.is_finite());
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("26 Cheyette/CMO/QRH/BaseCorr/EVaR")
        );
        ws.set_scenario(ScenarioKind::RatesPlus50Bp);
        let stressed = ws.step_tick();
        assert_ne!(
            stressed.cheyette_cmo_qrh_evar_metric,
            frame.cheyette_cmo_qrh_evar_metric
        );
    }
}

#[cfg(test)]
mod domain_suite_29_tests {
    use super::*;

    #[test]
    fn suite29_tab_metric_dashboard_and_rates_reactivity() {
        assert_eq!(TabId::ALL.len(), 77);
        assert!(TabId::ALL.contains(&TabId::AswMf1fAndreasenHugeBessEp));
        assert_eq!(
            TabId::AswMf1fAndreasenHugeBessEp.label(),
            "27 ASW / MF1F / AH-LV / BESS / Idzorek-DP"
        );
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert!(frame.asw_mf1f_ah_bess_ep_metric.is_finite());
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("27 ASW/MF1F/AH-LV/BESS/Idzorek-DP")
        );
        let normal_next = ws.step_tick();
        assert_ne!(
            normal_next.asw_mf1f_ah_bess_ep_metric,
            frame.asw_mf1f_ah_bess_ep_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.asw_mf1f_ah_bess_ep_metric,
                normal_next.asw_mf1f_ah_bess_ep_metric
            );
        }
    }
}

#[cfg(test)]
mod domain_suite_30_tests {
    use super::*;
    #[test]
    fn suite30_tab_metric_dashboard_and_stress_reactivity() {
        assert_eq!(TabId::ALL.len(), 77);
        assert_eq!(
            TabId::MtmXccsColvaZabrPdvWeatherEnb.label(),
            "28 MtM-XCCS / CTD-COLVA / ZABR-PDV / Weather / ENB"
        );
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        assert!(normal.mtmxccs_colva_zabr_weather_enb_metric.is_finite());
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("28 MtM-XCCS/COLVA/ZABR/Weather/ENB")
        );
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("BL / Risk Parity ERC")
        );
        let next = ws.step_tick();
        assert_ne!(
            next.mtmxccs_colva_zabr_weather_enb_metric,
            normal.mtmxccs_colva_zabr_weather_enb_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.mtmxccs_colva_zabr_weather_enb_metric,
                normal.mtmxccs_colva_zabr_weather_enb_metric
            );
        }
    }
}

#[cfg(test)]
mod domain_suite_31_tests {
    use super::*;

    #[test]
    fn suite31_tab29_terminal_renders_and_ticks_cleanly() {
        assert_eq!(TabId::ALL.len(), 77);
        assert_eq!(
            TabId::BurnoutCmcdsQuinticAfvCsSchur.label(),
            "29 Burnout-MBS / CMCDS / Quintic-OU / AFV-CB / CS-Swaption / Schur-HRP"
        );
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        assert!(normal.burnout_cmcds_quintic_afv_cs_schur_metric.is_finite());
        let dashboard = ws.render_dashboard(&normal, false);
        assert!(dashboard.contains("29 Burnout/CMCDS/Quintic/AFV/Schur"));
        assert!(dashboard.contains("BL / Risk Parity ERC"));
        let next = ws.step_tick();
        assert_ne!(
            next.burnout_cmcds_quintic_afv_cs_schur_metric,
            normal.burnout_cmcds_quintic_afv_cs_schur_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.burnout_cmcds_quintic_afv_cs_schur_metric,
                normal.burnout_cmcds_quintic_afv_cs_schur_metric
            );
        }
    }
}

#[cfg(test)]
mod domain_suite_32_tests {
    use super::*;

    #[test]
    fn suite32_tab_metric_dashboard_and_stress_reactivity() {
        assert_eq!(TabId::ALL.len(), 77);
        assert!(TabId::ALL.contains(&TabId::SinHouMbRBergomiChoiAfsExpectile));
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        assert!(
            normal
                .sinhou_mb_rbergomi_choi_afs_expectile_metric
                .is_finite()
        );
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("30 SinhOU/MB/mrBergomi/Choi/EXVaR")
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert!(
                stressed
                    .sinhou_mb_rbergomi_choi_afs_expectile_metric
                    .is_finite()
            );
            assert_ne!(
                stressed.sinhou_mb_rbergomi_choi_afs_expectile_metric,
                normal.sinhou_mb_rbergomi_choi_afs_expectile_metric
            );
        }
    }
}

#[cfg(test)]
mod domain_suite_33_tests {
    use super::*;

    #[test]
    fn suite33_terminal_reactivity() {
        assert_eq!(TabId::ALL.len(), 77);
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        let next = ws.step_tick();
        assert_ne!(
            next.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric,
            normal.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric,
                normal.monotone_convex_clark_saccr_rvol_cartea_paleologo_metric
            );
        }
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("31 HW-MC/Clark/SACCR/rVol/Cartea/Paleo")
        );
        assert!(
            ws.render_dashboard(&normal, false)
                .contains("BL / Risk Parity ERC")
        );
    }
}

#[cfg(test)]
mod domain_suite_34_tests {
    use super::*;
    #[test]
    fn suite34_tab32_metric_renders() {
        let mut ws = DeskWorkstation::new();
        let frame = ws.step_tick();
        assert_eq!(
            TabId::KrdLpiCocoBergomi2fCgmyGueantMeucci.label(),
            "32 KRD / LPI / CoCo / Bergomi / CGMY / AAD"
        );
        assert!(
            frame
                .krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
                .is_finite()
        );
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("32 KRD/LPI/CoCo/Bergomi2F/CGMY/Guéant/Meucci")
        );
        let next = ws.step_tick();
        assert_ne!(
            next.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric,
            frame.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
        );
        for scenario in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let stressed = ws.step_tick();
            assert_ne!(
                stressed.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric,
                frame.krd_lpi_coco_bergomi2f_cgmy_gueant_meucci_metric
            );
        }
    }
}
#[cfg(test)]
mod domain_suite_35_tests {
    use super::*;
    #[test]
    fn suite35_tab33_metric_renders_and_reacts() {
        let mut ws = DeskWorkstation::new();
        let f = ws.step_tick();
        assert_eq!(
            TabId::RkhsFmmVvLsvJltCtdQrhParisianAidNco.label(),
            "33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML"
        );
        assert!(
            ws.render_dashboard(&f, false)
                .contains("33 RKHS/FMM/VV/LSV/JLT/CTD/QRH/NCO/DML live")
        );
        for sc in [
            ScenarioKind::RatesPlus50Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::CreditWideningPlus100Bp,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(sc);
            let g = ws.step_tick();
            assert!(
                g.rkhs_fmm_vv_lsv_jlt_ctd_qrh_parisian_aid_propagator_nco_dml_metric
                    .is_finite()
            );
        }
    }
}
#[cfg(test)]
mod domain_suite_36_tests {
    use super::*;
    #[test]
    fn suite36_tab_metric_tick_and_rate_reactivity() {
        assert_eq!(TabId::ALL.len(), 77);
        assert!(
            TabId::MbsIoPoRmvMoXvaDilutionSigGarchPot
                .label()
                .starts_with("34 MBS-IO-PO")
        );
        let mut ws = DeskWorkstation::new();
        let a = ws.step_tick();
        assert!(a.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric.is_finite());
        ws.set_scenario(ScenarioKind::RatesPlus50Bp);
        let b = ws.step_tick();
        assert!(b.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric.is_finite());
        assert_ne!(
            a.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric,
            b.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric
        );
    }
    #[test]
    fn suite36_tab34_terminal_renders_and_reacts() {
        let mut ws = DeskWorkstation::new();
        assert_eq!(
            TabId::MbsIoPoRmvMoXvaDilutionSigGarchPot.label(),
            "34 MBS-IO-PO / RMV / MO-XVA / Dilution / Sig / GARCH / POT / RIE / PIN / GP / Brinson / IERP"
        );
        let mut vals = Vec::new();
        for scenario in [
            ScenarioKind::RatesParallelPlus50Bp,
            ScenarioKind::VolShockPlus5Pts,
            ScenarioKind::CreditSpreadWidening100Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let f = ws.step_tick();
            vals.push(f.mbs_rmv_mo_sig_garch_rie_gp_ierp_metric);
            let frame = ws.render_dashboard(&f, false);
            assert!(frame.contains("34 MBS/RMV/MO/CB/Sig/HN/POT/RIE/PIN/GP/BC/IERP"));
            assert!(!frame.contains('\u{1b}') || frame.contains("34"));
        }
        assert!(vals.iter().all(|x| x.is_finite()));
        assert!(vals.windows(2).any(|w| (w[1] - w[0]).abs() > 1e-12));
    }
}

#[cfg(test)]
mod domain_suite_37_tests {
    use super::*;
    #[test]
    fn suite37_tab35_metric_reacts_across_scenarios() {
        let mut ws = DeskWorkstation::new();
        let mut vals = Vec::new();
        for scenario in [
            ScenarioKind::RatesParallelPlus50Bp,
            ScenarioKind::VolShockPlus5Pts,
            ScenarioKind::CreditSpreadWidening100Bp,
            ScenarioKind::EquityCrashMinus15Pct,
            ScenarioKind::StagflationPlus150Bp,
        ] {
            ws.set_scenario(scenario);
            let frame = ws.step_tick();
            assert!(
                frame
                    .kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric
                    .is_finite()
            );
            assert!(ws.render_dashboard(&frame, false).contains(
                "35 Kalman/VECM/HJ/IFRS9/CPT/Hasbrouck/FV/Lando/LPPL/CB-Arb/Causal/Malliavin"
            ));
            vals.push(
                frame
                    .kalman_vecm_hj_ifrs9_cpt_hasbrouck_fv_lando_lppl_cbarb_causal_malliavin_metric,
            );
        }
        assert!(vals.windows(2).any(|w| (w[1] - w[0]).abs() > 1e-12));
    }
}

#[cfg(test)]
mod domain_suite_38_tests {
    use super::*;
    #[test]
    fn suite38_terminal_engine_dual_fd_parity() {
        use ql_ad::Dual;
        fn make<T: ql_ad::Real>(x: T) -> BergomiForwardStartNapoleonVolSpec<T> {
            BergomiForwardStartNapoleonVolSpec {
                initial_forward_variance: x * x,
                vol_of_vol_omega: T::from_f64(0.3),
                kappa_short_k1: T::from_f64(3.0),
                kappa_long_k2: T::from_f64(0.5),
                long_factor_weight_theta: T::from_f64(0.5),
                factor_corr_rho12: T::from_f64(0.2),
                spot_vol_corr_rho_s1: T::from_f64(-0.5),
                spot_vol_corr_rho_s2: T::from_f64(-0.3),
                forward_start_t1: T::from_f64(0.5),
                forward_end_t2: T::one(),
                forward_moneyness_k: T::one(),
                napoleon_base_coupon: T::from_f64(0.01),
                napoleon_global_floor: T::zero(),
                napoleon_subperiods_m: 4,
            }
        }
        let h = 1e-5;
        let fd = (evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&make(0.2 + h))
            .forward_atm_implied_vol
            - evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&make(0.2 - h))
                .forward_atm_implied_vol)
            / (2.0 * h);
        let dual = make(Dual::<1>::variable(0.2, 0));
        let ad = evaluate_bergomi_forward_start_smile_and_napoleon_convexity(&dual)
            .forward_atm_implied_vol;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn suite38_tab36_terminal_renders_and_reacts() {
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        let normal_metric = normal.mf_rates_wildcard_ctd_pivot_tarf_bergomi_metric;
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let stressed = ws.step_tick();
        assert!(
            stressed
                .mf_rates_wildcard_ctd_pivot_tarf_bergomi_metric
                .is_finite()
        );
        assert!(
            (stressed.mf_rates_wildcard_ctd_pivot_tarf_bergomi_metric - normal_metric).abs()
                > 1e-12
        );
        assert!(ws.render_dashboard(&stressed, false).contains(
            "36 MF-Rates/Wildcard/TARF/Bergomi/Schönbucher/FRTB-CVA/Lookback/VTA/FLAM/ACD/HLZ/MLMC"
        ));
        assert_eq!(crate::TabId::TAB_COUNT, 77);
        assert_eq!(TabId::ALL.len(), 77);
        assert!(TabId::ALL.contains(&TabId::MfRatesWildcardCtdPivotTarfBergomi));
    }
}

#[cfg(test)]
mod domain_suite_39_tests {
    use super::*;

    #[test]
    fn suite39_tab37_terminal_renders_and_reacts() {
        let mut ws = DeskWorkstation::new();
        let normal = ws.step_tick();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let stressed = ws.step_tick();
        assert!(
            stressed
                .cochrane_ftpl_g2_twist_essvi_offshore_cb_metric
                .is_finite()
        );
        assert!(
            (stressed.cochrane_ftpl_g2_twist_essvi_offshore_cb_metric
                - normal.cochrane_ftpl_g2_twist_essvi_offshore_cb_metric)
                .abs()
                > 1e-12
        );
        assert!(
            ws.render_dashboard(&stressed, false).contains(
                "37 FTPL/G2/eSSVI/CB/CoCo/MPOR/Double-Heston/Curran/Warin/Roll/Busseti/KKT"
            )
        );
        assert_eq!(crate::TabId::TAB_COUNT, 77);
        assert_eq!(TabId::ALL.len(), 77);
        assert!(TabId::ALL.contains(&TabId::CochraneFtplG2TwistEssviOffshoreCb));
    }

    #[test]
    fn suite39_terminal_dual_fd_parity() {
        use ql_ad::Dual;
        fn eval<T: ql_ad::Real>(x: T) -> T {
            ql_inflation::evaluate_cochrane_ftpl_and_maturity_inflation_smoothing(
                &ql_inflation::CochraneFiscalTheoryPriceLevelSpec {
                    debt_maturities_years: vec![T::one()],
                    nominal_debt_face_values: vec![x],
                    nominal_zero_yields: vec![T::from_f64(0.03)],
                    real_zero_yields: vec![T::from_f64(0.01)],
                    primary_real_surpluses: vec![T::from_f64(5.0)],
                    terminal_surplus_growth_rate: T::from_f64(0.02),
                    fiscal_surplus_shock_fraction: T::from_f64(-0.05),
                    expected_future_inflation_rise: T::from_f64(0.01),
                    inflation_risk_beta: T::from_f64(0.4),
                    inflation_volatility: T::from_f64(0.1),
                    tips_liquidity_spread: T::from_f64(0.002),
                },
            )
            .equilibrium_price_level_p0
        }
        let h = 1e-5;
        let fd = (eval(100.0 + h) - eval(100.0 - h)) / (2.0 * h);
        let d = eval(Dual::<1>::variable(100.0, 0));
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_40_tests {
    use super::*;
    #[test]
    fn suite40_tab38_terminal_renders_and_reacts() {
        let mut ws = DeskWorkstation::default();
        let normal_frame = ws.step_tick();
        let normal = ws.render_dashboard(&normal_frame, false);
        ws.scenario = ScenarioKind::CreditWideningPlus100Bp;
        let stressed_frame = ws.step_tick();
        let stressed = ws.render_dashboard(&stressed_frame, false);
        assert_eq!(TabId::ALL.len(), 77);
        assert!(TabId::ALL.contains(&TabId::BcbsIrrbbTurfusRfrBaranMtmBourgeyBg));
        assert!(
            normal
                .contains("38 IRRBB/RFR/XCCY/BG/CoCo/SIMM/rBergomi/Spread/BESS/Dark/Risk/Diff-PCA")
        );
        assert!(
            stressed
                .contains("38 IRRBB/RFR/XCCY/BG/CoCo/SIMM/rBergomi/Spread/BESS/Dark/Risk/Diff-PCA")
        );
        assert!(
            (stressed_frame.bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_metric
                - normal_frame.bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_metric)
                .abs()
                > 1e-12
        );
        assert!(
            (ws.bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_live_metric()
                - DeskWorkstation::bcbs_irrbb_turfus_rfr_baran_mtm_bourgey_bg_live_metric_from_market_with_bump(
                    &ws.scenario_market(),
                    0.01,
                ))
            .abs()
                < 1e-12
        );
    }
}

#[cfg(test)]
mod domain_suite_41_tests {
    use super::*;
    #[test]
    fn suite41_tab39_terminal_renders_and_reacts() {
        assert_eq!(crate::TabId::TAB_COUNT, 77);
        assert!(TabId::VeronesiConvLandoGenGatheralCorridor.label().starts_with("39 Veronesi / Lando-Log / Corridor-Var / Russian / Bates-SVJ / Koziol-CB / MFG / BCBS-Floor / CJP-LOB / ETF-AP / IC-Decay / Glasserman-IS"));
        let mut ws = DeskWorkstation::default();
        let a = ws.step_tick();
        let b = ws.step_tick();
        assert!(
            a.veronesi_conv_lando_gen_gatheral_corridor_metric
                .is_finite()
                && b.veronesi_conv_lando_gen_gatheral_corridor_metric
                    .is_finite()
        );
        assert!(a.veronesi_conv_lando_gen_gatheral_corridor_metric != 0.0);
        assert!(
            ws.render_dashboard(&a, false)
                .contains("39 Veronesi/Lando/Corridor/Russian/Bates/Koziol/MFG/BCBS/CJP/ETF/IC/IS")
        );
    }
}

#[cfg(test)]
mod domain_suite_42_tests {
    use super::*;

    #[test]
    fn suite42_tab40_terminal_renders_and_reacts() {
        assert_eq!(crate::TabId::TAB_COUNT, 77);
        assert_eq!(TabId::ALL.len(), 77);
        assert!(
            TabId::RomerTaylorDamodaranCrpChristoffersenVar
                .label()
                .starts_with("40 Taylor-Curve / Damodaran-CRP / Cumulant-VaR / LOB-LVaR / Prospect-CGO / Regret-Match / Beneish-HKAS32 / Campisi-FI / Heckman-GMM / Elastic-Net / Reset-CB / FSB-BailIn")
        );
        let mut ws = DeskWorkstation::new();
        ws.set_scenario(ScenarioKind::Normal);
        let normal = ws.step_tick();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let stressed = ws.step_tick();
        assert!(
            normal
                .romer_taylor_damodaran_crp_christoffersen_var_metric
                .is_finite()
        );
        assert!(
            stressed
                .romer_taylor_damodaran_crp_christoffersen_var_metric
                .is_finite()
        );
        assert_ne!(
            normal.romer_taylor_damodaran_crp_christoffersen_var_metric,
            stressed.romer_taylor_damodaran_crp_christoffersen_var_metric
        );
        assert!(ws.render_dashboard(&stressed, false).contains(
            "40 Taylor/Damodaran/VaR/LVaR/CGO/Regret/Beneish/Campisi/Heckman/Elastic/ResetCB/FSB"
        ));
    }
}

#[cfg(test)]
mod domain_suite_43_tests {
    use super::*;

    #[test]
    fn suite43_tab41_terminal_renders_and_reacts() {
        assert_eq!(crate::TabId::TAB_COUNT, 77);
        assert_eq!(TabId::ALL.len(), 77);
        assert!(
            TabId::FabozziKrdCmoPacAsquithCbAsc470Colva
                .label()
                .starts_with("41 KRD-LSC / CMO-PAC / Asquith-CB / ASU2024-CPDI / ASC470-Mod / ColVA-CTD / Hill-EVT / Bouchaud-TLF / Karnosky-FX / DGTW-Timing / Brock-Hommes / Penman-ReOI")
        );
        let mut ws = DeskWorkstation::new();
        ws.set_scenario(ScenarioKind::Normal);
        let normal = ws.step_tick();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let credit = ws.step_tick();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let stag = ws.step_tick();
        assert!(
            normal
                .fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric
                .is_finite()
        );
        assert!(
            (credit.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric
                - normal.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric)
                .abs()
                > 1e-12
        );
        assert!(
            (stag.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric
                - normal.fabozzi_krd_cmo_pac_asquith_cb_asc470_colva_metric)
                .abs()
                > 1e-12
        );
        assert!(
            ws.render_dashboard(&stag, false)
                .contains("41 KRD-LSC/CMO-PAC/Asquith-CB/ASU2024/ASC470/ColVA/Hill/Bouchaud/Karnosky/DGTW/Brock/Penman")
        );
    }
}

#[cfg(test)]
mod domain_suite_44_tests {
    use super::*;
    #[test]
    fn suite44_tab42_terminal_renders_and_reacts() {
        assert_eq!(crate::TabId::TAB_COUNT, 77);
        assert_eq!(TabId::ALL.len(), 77);
        assert!(
            TabId::CheyetteCapletAdbNdirsZabrQrhCmcdsAdi
                .label()
                .starts_with("42 Cheyette-LV / ADB-NDIRS / Andreasen-ZABR / QRH-Zumbach / CMCDS-Morini / Turlakov-WWR / HHW-ADI / CappedCall-MW / Ahmadi-EVaR / Schur-HRP / BEKK-GARCH / Schilit-Sloan")
        );
        let mut ws = DeskWorkstation::default();
        let m0 = ws.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_live_metric();
        let f1 = ws.step_tick();
        assert!((f1.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_metric - m0).abs() > 1e-12);
        ws.set_scenario(ScenarioKind::Normal);
        let normal = ws.step_tick();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let credit = ws.step_tick();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let stag = ws.step_tick();
        assert!(
            (credit.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_metric
                - normal.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_metric)
                .abs()
                > 1e-12
        );
        assert!(
            (stag.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_metric
                - normal.cheyette_caplet_adb_ndirs_zabr_qrh_cmcds_adi_metric)
                .abs()
                > 1e-12
        );
        assert!(
            ws.render_dashboard(&stag, false)
                .contains("42 Cheyette/NDIRS/ZABR/QRH/CMCDS/WWR/HHW/CB/EVaR/Schur/BEKK/Quality")
        );
    }
}

#[cfg(test)]
mod domain_suite_45_tests {
    use super::*;
    #[test]
    fn suite45_tab43_terminal_renders_and_reacts() {
        assert_eq!(crate::TabId::TAB_COUNT, 77);
        assert_eq!(TabId::ALL.len(), 77);
        assert!(
            TabId::RebonatoCmsLelandToftAscotPoolRollRie
                .label()
                .starts_with("43 Rebonato-CMS / Leland-Toft / Calamos-Ascot / SpecPool-Roll / Glosten-Milgrom / Potters-RIE / EGARCH-Hamilton / Wooldridge-2SLS / McNeil-Euler / Bacon-Omega / Damodaran-RD / Arrow-Debreu")
        );
        let mut ws = DeskWorkstation::new();
        let a = ws.rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric();
        let _ = ws.step_tick();
        assert!(
            ws.rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric()
                .is_finite()
        );
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.rebonato_cms_leland_toft_ascot_pool_roll_rie_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("43 Rebonato-CMS/Leland-Toft/Calamos/RIE/Arrow-Debreu")
        );
    }
}

#[cfg(test)]
mod domain_suite_46_tests {
    use super::*;
    #[test]
    fn suite46_tab44_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::TuckmanRepoBkFvaDecsIoPoKylePinGjrLppl
                .label()
                .starts_with("44 Tuckman-Repo")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.tuckman_repo_bk_fva_decs_io_po_kyle_pin_gjr_lppl_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("44 Tuckman-Repo")
        );
    }
}

#[cfg(test)]
mod domain_suite_47_tests {
    use super::*;
    #[test]
    fn suite47_tab45_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::HaganWestCtdSwitchAfvCbKthDefaultSaCva
                .label()
                .starts_with("45 Hagan-West-MC")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.hagan_west_ctd_switch_afv_cb_kth_default_sa_cva_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("45 Hagan-West-MC")
        );
    }
}

#[cfg(test)]
mod domain_suite_48_tests {
    use super::*;
    #[test]
    fn suite48_tab46_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::BrigoCsaCtdCallableOasJltCorridorLsm
                .label()
                .starts_with("46 Multi-CSA-CTD")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.brigo_csa_ctd_callable_oas_jlt_corridor_lsm_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("46 Multi-CSA-CTD")
        );
    }
}

#[cfg(test)]
mod domain_suite_49_tests {
    use super::*;
    #[test]
    fn suite49_tab47_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::BrigoNssKrdTipsFloorRmvTriangleStickyLv
                .label()
                .starts_with("47 NSS-KRD")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.brigo_nss_krd_tips_floor_rmv_triangle_sticky_lv_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(ws.render_dashboard(&frame, false).contains("47 NSS-KRD"));
    }
}

#[cfg(test)]
mod domain_suite_50_tests {
    use super::*;
    #[test]
    fn suite50_tab48_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::HaganWestTensionPacCmoFamaMacBethAsviBgk
                .label()
                .starts_with("48 Tension-Curve")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.hagan_west_tension_pac_cmo_fama_mac_beth_asvi_bgk_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("48 Tension-Curve")
        );
    }
}

#[cfg(test)]
mod domain_suite_51_tests {
    use super::*;
    #[test]
    fn suite51_tab49_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::BrigoXccyMtmSpecPoolMarshallOlkinVixVvix
                .label()
                .starts_with("49 XCCY-MtM")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.brigo_xccy_mtm_spec_pool_marshall_olkin_vix_vvix_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(ws.render_dashboard(&frame, false).contains("49 XCCY-MtM"));
    }
}

#[cfg(test)]
mod domain_suite_52_tests {
    use super::*;
    #[test]
    fn suite52_tab50_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::BrigoOisFallbackCollaredFrnCdsFepSparkPrdc
                .label()
                .starts_with("50 OIS-Fallback")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.brigo_ois_fallback_collared_frn_cds_fep_spark_prdc_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("50 OIS-Fallback")
        );
    }
}

#[cfg(test)]
mod domain_suite_53_tests {
    use super::*;
    #[test]
    fn suite53_tab51_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::AndersenCmsSpecialRepoCdsUpfrontSsrKvaCrack
                .label()
                .starts_with("51 CMS-Replication")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.andersen_cms_special_repo_cds_upfront_ssr_kva_crack_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(ws.render_dashboard(&frame, false).contains("51 CMS-Rep"));
    }
}

#[cfg(test)]
mod domain_suite_54_tests {
    use super::*;
    #[test]
    fn suite54_tab52_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::DarbyshireStepOisCbSoftCallCdxSkewMalzFly
                .label()
                .starts_with("52 Step-OIS-Turn")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.darbyshire_step_ois_cb_soft_call_cdx_skew_malz_fly_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("52 Step-OIS-Turn")
        );
    }
}

#[cfg(test)]
mod domain_suite_55_tests {
    use super::*;
    #[test]
    fn suite55_tab53_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::HenrardTenorBasisMsrRecaptureCdoBaseCorrKiko
                .label()
                .starts_with("53 Tenor-Basis-Mult")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.henrard_tenor_basis_msr_recapture_cdo_base_corr_kiko_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("53 Tenor-Basis/MSR-Recap")
        );
    }
}

#[cfg(test)]
mod domain_suite_56_tests {
    use super::*;
    #[test]
    fn suite56_tab54_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::BrigoDiffSwapCallableMuniBkmSkewBreakCva
                .label()
                .starts_with("54 Diff-Swap-Quanto")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.brigo_diff_swap_callable_muni_bkm_skew_break_cva_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("54 Diff-Swap/Muni-YTW")
        );
    }
}

#[cfg(test)]
mod domain_suite_57_tests {
    use super::*;
    #[test]
    fn suite57_tab55_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::AndersenSofrFutHayreMbsPedersenCdxDupireLv
                .label()
                .starts_with("55 SOFR-Fut-FRA")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.andersen_sofr_fut_hayre_mbs_pedersen_cdx_dupire_lv_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("55 SOFR-Fut/Hayre-MBS")
        );
    }
}

#[cfg(test)]
mod domain_suite_58_tests {
    use super::*;
    #[test]
    fn suite58_tab56_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::LittermanScheinkmanPcaFlyTfSplitCbDoubleTCdoLeeWings
                .label()
                .starts_with("56 PCA-Fly-Carry")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.litterman_scheinkman_pca_fly_tf_split_cb_double_t_cdo_lee_wings_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("56 PCA-Fly/TF-CB")
        );
    }
}

#[cfg(test)]
mod domain_suite_59_tests {
    use super::*;
    #[test]
    fn suite59_tab57_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::PiterbargCollateralOptTuckmanRolldownMertonKmvBartlettSabr
                .label()
                .starts_with("57 Piterbarg-Coll-Opt")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st =
            ws.piterbarg_collateral_opt_tuckman_rolldown_merton_kmv_bartlett_sabr_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("57 Piterbarg-Opt/Tuckman-Roll")
        );
    }
}

#[cfg(test)]
mod domain_suite_60_tests {
    use super::*;
    #[test]
    fn suite60_tab58_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::GswSvenssonCurveFabozziMbsEsfVasicekGordyIrbCarrWuVol
                .label()
                .starts_with("58 GSW-Svensson-Curve")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.gsw_svensson_curve_fabozzi_mbs_esf_vasicek_gordy_irb_carr_wu_vol_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("58 GSW-Svensson/MBS-ESF")
        );
    }
}

#[cfg(test)]
mod domain_suite_61_tests {
    use super::*;
    #[test]
    fn suite61_tab59_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::KimWrightShadowBrennanSchwartzCbDuffieSingletonRmvFukasawaGamma
                .label()
                .starts_with("59 Kim-Wright-Shadow")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws
            .kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric(
            );
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws
            .kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric(
            );
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws
            .kim_wright_shadow_brennan_schwartz_cb_duffie_singleton_rmv_fukasawa_gamma_live_metric(
            );
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("59 Kim-Wright-Shadow/Brennan-Schwartz-CB")
        );
    }
}

#[cfg(test)]
mod domain_suite_62_tests {
    use super::*;
    #[test]
    fn suite62_tab60_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::JpsUnspannedMacroKalotayRefundingHullWhiteNtdDermanSticky
                .label()
                .starts_with("60 JPS-Unspanned-Macro")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st =
            ws.jps_unspanned_macro_kalotay_refunding_hull_white_ntd_derman_sticky_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("60 JPS-Unspanned-Macro/Kalotay-Refunding")
        );
    }
}

#[cfg(test)]
mod domain_suite_63_tests {
    use super::*;
    #[test]
    fn suite63_tab61_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::LyashenkoFmmCapletTuckmanCpiSeasonalityBrigoAlfonsiCirFouqueSv
                .label()
                .starts_with("61 Lyashenko-FMM-Caplet")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws
            .lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws
            .lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws
            .lyashenko_fmm_caplet_tuckman_cpi_seasonality_brigo_alfonsi_cir_fouque_sv_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("61 Lyashenko-FMM-Caplet/Tuckman-CPI-Seasonality")
        );
    }
}

#[cfg(test)]
mod domain_suite_64_tests {
    use super::*;
    #[test]
    fn suite64_tab62_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::JszCanonicalGdtsmTsiveriotisFernandesCbLaurentGregoryCdoBroadieJainVarswap
                .label()
                .starts_with("62 JSZ-Canonical-GDTSM")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.jsz_canonical_gdtsm_tsiveriotis_fernandes_cb_laurent_gregory_cdo_broadie_jain_varswap_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("62 JSZ-Canonical-GDTSM/Tsiveriotis-Fernandes-CB")
        );
    }
}

#[cfg(test)]
mod domain_suite_65_tests {
    use super::*;
    #[test]
    fn suite65_tab63_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::BianchettiTwoCurveQuantoBdtCallableBondSchonbucherCdsSwaptionBatesSvj
                .label()
                .starts_with("63 Bianchetti-Two-Curve-Quanto")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.bianchetti_two_curve_quanto_bdt_callable_bond_schonbucher_cds_swaption_bates_svj_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("63 Bianchetti-Two-Curve-Quanto/BDT-Callable")
        );
    }
}

#[cfg(test)]
mod domain_suite_66_tests {
    use super::*;
    #[test]
    fn suite66_tab64_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::AndersenPiterbargCmsSpreadFabozziZBondCmoDuffieGarleanuAffineCdoCarrLee
                .label()
                .starts_with("64 Andersen-Piterbarg-CMS-Spread")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.andersen_piterbarg_cms_spread_fabozzi_z_bond_cmo_duffie_garleanu_affine_cdo_carr_lee_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("64 Andersen-Piterbarg-CMS-Spread/Fabozzi-Z-Bond")
        );
    }
}

#[cfg(test)]
mod domain_suite_67_tests {
    use super::*;
    #[test]
    fn suite67_tab65_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::PiterbargCtdCollateralTuckmanMakeWholeBondBrigoMoriniCmcdsHaganShiftedSabr
                .label()
                .starts_with("65 Piterbarg-CTD-Collateral")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.piterbarg_ctd_collateral_tuckman_make_whole_bond_brigo_morini_cmcds_hagan_shifted_sabr_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("65 Piterbarg-CTD-Collateral/Tuckman-Make-Whole")
        );
    }
}

#[cfg(test)]
mod domain_suite_68_tests {
    use super::*;
    #[test]
    fn suite68_tab66_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::MercurioFvaMultiCurveSwapIngersollContingentCallCbHullWhiteImpliedCopulaCdoGatheralJacquierSviJw
                .label()
                .starts_with("66 Mercurio-FVA-Swap")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.mercurio_fva_multi_curve_swap_ingersoll_contingent_call_cb_hull_white_implied_copula_cdo_gatheral_jacquier_svi_jw_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("66 Mercurio-FVA-Swap/Ingersoll-CB")
        );
    }
}

#[cfg(test)]
mod domain_suite_69_tests {
    use super::*;
    #[test]
    fn suite69_tab67_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::AmetranoBianchettiTenorBasisFabozziHayreIoPoStripHullBrigoWwrCdsBossensVannaVolgaBarrier
                .label()
                .starts_with("67 Ametrano-Bianchetti-Tenor-Basis")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.ametrano_bianchetti_tenor_basis_fabozzi_hayre_io_po_strip_hull_brigo_wwr_cds_bossens_vanna_volga_barrier_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("67 Ametrano-Bianchetti-Tenor-Basis/Fabozzi-Hayre-IO-PO")
        );
    }
}

#[cfg(test)]
mod domain_suite_70_tests {
    use super::*;
    #[test]
    fn suite70_tab68_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::SmithWilsonUfrKalotayWilliamsFabozziCallableSchonbucherClaytonCdoBuehlerVarianceCurve
                .label()
                .starts_with("68 Smith-Wilson-UFR")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.smith_wilson_ufr_kalotay_williams_fabozzi_callable_schonbucher_clayton_cdo_buehler_variance_curve_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("68 Smith-Wilson-UFR/Kalotay-Williams-Fabozzi")
        );
    }
}

#[cfg(test)]
mod domain_suite_71_tests {
    use super::*;
    #[test]
    fn suite71_tab69_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::McCullochTaxCurveAndersenBuffumCbOkaneRecoveryLockCdsDermanKaniImpliedTree
                .label()
                .starts_with("69 McCulloch-Tax-Curve")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.mcculloch_tax_curve_andersen_buffum_cb_okane_recovery_lock_cds_derman_kani_implied_tree_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("69 McCulloch-Tax-Curve/Andersen-Buffum-CB")
        );
    }
}

#[cfg(test)]
mod domain_suite_72_tests {
    use super::*;
    #[test]
    fn suite72_tab70_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::HaganWestMonotoneConvexFabozziEffectiveDurationOasFryeAltmanDownturnLgdGatheralJacquierSviArbitrage
                .label()
                .starts_with("70 Hagan-West-Monotone-Convex")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.hagan_west_monotone_convex_fabozzi_effective_duration_oas_frye_altman_downturn_lgd_gatheral_jacquier_svi_arbitrage_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("70 Hagan-West-Monotone-Convex/Fabozzi-Effective-Duration-OAS")
        );
    }
}

#[cfg(test)]
mod domain_suite_73_tests {
    use super::*;
    #[test]
    fn suite73_tab71_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::MercurioRussoIborOisBasisJarrowYildirimYoyInflationDuffieLandoHybridCdsBayerFrizGatheralRoughSkew
                .label()
                .starts_with("71 Mercurio-Russo-IBOR-OIS-Basis")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.mercurio_russo_ibor_ois_basis_jarrow_yildirim_yoy_inflation_duffie_lando_hybrid_cds_bayer_friz_gatheral_rough_skew_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("71 Mercurio-Russo-IBOR-OIS-Basis/Jarrow-Yildirim-YoY-Inflation")
        );
    }
}

#[cfg(test)]
mod domain_suite_74_tests {
    use super::*;
    #[test]
    fn suite74_tab72_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::FilipovicLarssonTrolleLinearRationalDuffieKrishnamurthySpecialRepoDuffieSaitaWangStochasticDefaultCarrMadanNeubergerLogContract
                .label()
                .starts_with("72 Filipovic-Larsson-Trolle-LRTSM")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.filipovic_larsson_trolle_linear_rational_duffie_krishnamurthy_special_repo_duffie_saita_wang_stochastic_default_carr_madan_neuberger_log_contract_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("72 Filipovic-Larsson-Trolle-LRTSM/Duffie-Krishnamurthy-Repo")
        );
    }
}

#[cfg(test)]
mod domain_suite_75_tests {
    use super::*;
    #[test]
    fn suite75_tab73_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::LyashenkoMercurioRfrInArrearsPennacchiSundaresanCocoAt1PedersenBrigoCdsIndexFepBroadieJainSeppSvjVolswap
                .label()
                .starts_with("73 Lyashenko-Mercurio-RFR")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.lyashenko_mercurio_rfr_in_arrears_pennacchi_sundaresan_coco_at1_pedersen_brigo_cds_index_fep_broadie_jain_sepp_svj_volswap_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("73 Lyashenko-Mercurio-RFR/Pennacchi-Sundaresan-CoCo")
        );
    }
}

#[cfg(test)]
mod domain_suite_76_tests {
    use super::*;
    #[test]
    fn suite76_tab74_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::CieslakPovalaTrendInflationCycleGabaixKrishnamurthyVigneronMbsOasErraisGieseckeGoldbergHawkesCreditMartinBondarenkoSvixLowerBound
                .label()
                .starts_with("74 Cieslak-Povala-Cycle")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.cieslak_povala_trend_inflation_cycle_gabaix_krishnamurthy_vigneron_mbs_oas_errais_giesecke_goldberg_hawkes_credit_martin_bondarenko_svix_lower_bound_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("74 Cieslak-Povala-Cycle/Gabaix-Krishnamurthy-Vigneron-MBS")
        );
    }
}

#[cfg(test)]
mod domain_suite_77_tests {
    use super::*;
    #[test]
    fn suite77_tab75_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::BianchettiMercurioMultiCurveFraQuantoAcharyaAmihudBharathLiquidityStressBrigoAlfonsiCorrelatedCirCdsFukasawaAsymptoticSkewVarswap
                .label()
                .starts_with("75 Bianchetti-Mercurio-FRA-Quanto")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.bianchetti_mercurio_multi_curve_fra_quanto_acharya_amihud_bharath_liquidity_stress_brigo_alfonsi_correlated_cir_cds_fukasawa_asymptotic_skew_varswap_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("75 Bianchetti-Mercurio-FRA-Quanto/Acharya-Amihud-Bharath-Liq")
        );
    }
}

#[cfg(test)]
mod domain_suite_78_tests {
    use super::*;
    #[test]
    fn suite78_tab76_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::HansonSteinGreenwoodVayanosDurationSupplyFleckensteinLongstaffLustigTipsPuzzleCollinDufresneGoldsteinYangSuperSeniorCdoTodorovOptionSpotVolatility
                .label()
                .starts_with("76 Hanson-Stein-Duration")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.hanson_stein_greenwood_vayanos_duration_supply_fleckenstein_longstaff_lustig_tips_puzzle_collin_dufresne_goldstein_yang_super_senior_cdo_todorov_option_spot_volatility_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("76 Hanson-Stein-Duration/Fleckenstein-Longstaff-Lustig-TIPS")
        );
    }
}

#[cfg(test)]
mod domain_suite_79_tests {
    use super::*;
    #[test]
    fn suite79_tab77_terminal_renders_and_reacts() {
        assert_eq!(TabId::TAB_COUNT, 77);
        assert!(
            TabId::KlinglerSyrstadSofrRepoSpikeLongstaffMithalNeisCdsBondBasisPanSingletonSovereignCdsMedvedevScailletShortMaturitySmile
                .label()
                .starts_with("77 Klingler-Syrstad-SOFR-Spike")
        );
        let mut ws = DeskWorkstation::default();
        let a = ws.klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric();
        ws.set_scenario(ScenarioKind::CreditWideningPlus100Bp);
        let c = ws.klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric();
        ws.set_scenario(ScenarioKind::StagflationPlus150Bp);
        let st = ws.klingler_syrstad_sofr_repo_spike_longstaff_mithal_neis_cds_bond_basis_pan_singleton_sovereign_cds_medvedev_scaillet_short_maturity_smile_live_metric();
        assert!(a.is_finite() && (a - c).abs() > 1e-12 && (c - st).abs() > 1e-12);
        let frame = ws.step_tick();
        assert!(
            ws.render_dashboard(&frame, false)
                .contains("77 Klingler-Syrstad-SOFR-Spike/Longstaff-Mithal-Neis-CDS-Bond")
        );
    }
}
