#pragma once



#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace solution_audit_output
{
	inline bool needs_header(const std::string& path)
	{
		std::ifstream in(path, std::ios::binary);
		return !in.good() || in.peek() == std::char_traits<char>::eof();
	}

	// Append one aggregate row and one row per route. Use fresh output files (or
	// remove old files) when changing the column layout so rows remain consistent.
	inline void append_solution_stats(
		const std::string& summary_path,
		const std::string& routes_path,
		const std::string& instance_name,
		int repetition,
		int bestknown,
		int score,
		double runtime_seconds,
		double gap_percent,
		const SolutionPerfStats& audit)
	{
		const bool summary_header = needs_header(summary_path);
		std::ofstream summary(summary_path, std::ios::out | std::ios::app);
		if (!summary)
			throw std::runtime_error("Could not open solution audit file: " + summary_path);

		if (summary_header)
		{
			summary
				<< "instance;rep;bestknown;score;runtime_s;gap_pct;"
				<< "served_customers;nonempty_routes;route_count;"
				<< "avg_route_duration_h;avg_route_time_util_pct;max_route_time_util_pct;"
				<< "routes_time_util_ge_90pct;pct_routes_time_util_ge_90pct;"
				<< "avg_weight_util_pct;max_weight_util_pct;routes_weight_util_ge_90pct;"
				<< "pct_routes_weight_util_ge_90pct;"
				<< "avg_volume_util_pct;max_volume_util_pct;routes_volume_util_ge_90pct;"
				<< "pct_routes_volume_util_ge_90pct;"
				<< "min_customer_tw_slack_h;avg_customer_tw_slack_h;min_max_shift_h;"
				<< "avg_waiting_time_customer_h;pct_break_at_end_depot;"
				<< "avg_break_start_after_EDT_h;avg_break_position_norm;avg_score_per_customer\n";
		}

		summary << std::setprecision(12)
			<< instance_name << ';' << repetition << ';' << bestknown << ';' << score << ';'
			<< runtime_seconds << ';' << gap_percent << ';'
			<< audit.total_served_customers << ';'
			<< audit.nonempty_route_count << ';' << audit.route_count << ';'
			<< audit.avg_route_duration << ';'
			<< 100.0 * audit.avg_route_time_utilization << ';'
			<< 100.0 * audit.max_route_time_utilization << ';'
			<< audit.routes_time_util_ge_90pct << ';'
			<< audit.pct_routes_time_util_ge_90pct << ';'
			<< 100.0 * audit.avg_weight_utilization << ';'
			<< 100.0 * audit.max_weight_utilization << ';'
			<< audit.routes_weight_util_ge_90pct << ';'
			<< audit.pct_routes_weight_util_ge_90pct << ';'
			<< 100.0 * audit.avg_volume_utilization << ';'
			<< 100.0 * audit.max_volume_utilization << ';'
			<< audit.routes_volume_util_ge_90pct << ';'
			<< audit.pct_routes_volume_util_ge_90pct << ';'
			<< audit.min_customer_tw_slack << ';'
			<< audit.avg_customer_tw_slack << ';'
			<< audit.min_max_shift << ';'
			<< audit.avg_waiting_time_per_customer << ';'
			<< audit.pct_break_at_end_depot << ';'
			<< audit.avg_break_start_time << ';'
			<< audit.avg_break_position_norm << ';'
			<< audit.avg_score_per_served_customer << '\n';
		summary.close();

		const bool routes_header = needs_header(routes_path);
		std::ofstream routes(routes_path, std::ios::out | std::ios::app);
		if (!routes)
			throw std::runtime_error("Could not open route audit file: " + routes_path);

		if (routes_header)
		{
			routes
				<< "instance;rep;route;score;served_customers;duration_h;time_budget_h;"
				<< "time_util_pct;weight_used;weight_capacity;weight_util_pct;weight_slack;"
				<< "volume_used;volume_capacity;volume_util_pct;volume_slack;"
				<< "min_customer_tw_slack_h;avg_customer_tw_slack_h;"
				<< "customers_tw_slack_le_0_1h;min_max_shift_h;"
				<< "break_taken;break_at_end_depot;break_start_after_EDT_h;"
				<< "break_start_slack_h;break_position_norm;"
				<< "total_waiting_h;avg_waiting_per_customer_h\n";
		}

		routes << std::setprecision(12);
		for (const RoutePerfStats& route : audit.route_stats)
		{
			routes
				<< instance_name << ';' << repetition << ';' << route.route_index << ';'
				<< route.collected_score << ';' << route.served_customers << ';'
				<< route.route_duration << ';' << route.route_time_budget << ';'
				<< 100.0 * route.route_time_utilization << ';'
				<< route.weight_used << ';' << route.weight_capacity << ';'
				<< 100.0 * route.weight_utilization << ';' << route.weight_slack << ';'
				<< route.volume_used << ';' << route.volume_capacity << ';'
				<< 100.0 * route.volume_utilization << ';' << route.volume_slack << ';'
				<< route.min_customer_tw_slack << ';'
				<< route.avg_customer_tw_slack << ';'
				<< route.customers_tw_slack_le_0_1 << ';'
				<< route.min_max_shift << ';'
				<< (route.break_taken ? 1 : 0) << ';'
				<< (route.break_at_end_depot ? 1 : 0) << ';'
				<< route.break_start_time << ';' << route.break_start_slack << ';'
				<< route.break_position_norm << ';' << route.total_waiting_time << ';'
				<< route.avg_waiting_time_per_customer << '\n';
		}
	}
}

