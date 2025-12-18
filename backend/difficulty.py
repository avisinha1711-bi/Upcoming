# backend/adaptive_difficulty.py
class AdaptiveDifficultySystem:
    """AI that adjusts game difficulty based on player performance"""
    
    def __init__(self):
        self.player_history = {}
        self.difficulty_profiles = {
            'beginner': {'speed_interval': 10, 'speed_increase': 0.2, 'max_bubbles': 8},
            'casual': {'speed_interval': 7, 'speed_increase': 0.3, 'max_bubbles': 12},
            'pro': {'speed_interval': 5, 'speed_increase': 0.5, 'max_bubbles': 15},
            'expert': {'speed_interval': 3, 'speed_increase': 0.7, 'max_bubbles': 20}
        }
    
    def analyze_player(self, player_id, game_data):
        """Analyze player performance and adjust difficulty"""
        
        accuracy = game_data['hits'] / max(game_data['shots'], 1)
        avg_survival_time = game_data['total_time'] / max(game_data['games_played'], 1)
        max_speed_level = game_data['max_speed_level']
        
        # Machine Learning-based classification
        if accuracy > 0.85 and max_speed_level > 10:
            profile = 'expert'
        elif accuracy > 0.70 and max_speed_level > 6:
            profile = 'pro'
        elif accuracy > 0.50:
            profile = 'casual'
        else:
            profile = 'beginner'
        
        # Adjust for frustration detection
        if game_data.get('rage_quits', 0) > 2:
            profile = self._ease_difficulty(profile)
        
        return self.difficulty_profiles[profile]
    
    def generate_personalized_challenge(self, player_id):
        """Create custom bubble patterns for player"""
        weaknesses = self._identify_weaknesses(player_id)
        
        if 'left_side' in weaknesses:
            return self._create_left_heavy_pattern()
        elif 'fast_bubbles' in weaknesses:
            return self._create_speed_challenge_pattern()
        
        return self._create_random_pattern()
