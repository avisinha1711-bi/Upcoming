# backend/content_generator.py
import random
from PIL import Image, ImageDraw
import io
import base64

class DynamicContentGenerator:
    """Generate game content dynamically"""
    
    def generate_bubble_themes(self, player_level):
        """Create themed bubble sets based on player progress"""
        themes = {
            1: {'name': 'Classic', 'colors': ['#FF5252', '#FF4081', '#E040FB']},
            5: {'name': 'Ocean', 'colors': ['#2196F3', '#03A9F4', '#00BCD4']},
            10: {'name': 'Forest', 'colors': ['#4CAF50', '#8BC34A', '#CDDC39']},
            20: {'name': 'Galaxy', 'colors': ['#9C27B0', '#673AB7', '#3F51B5']}
        }
        
        # Unlock themes as player progresses
        unlocked = [themes[lvl] for lvl in themes if lvl <= player_level]
        return unlocked
    
    def create_special_bubbles(self, special_event):
        """Generate special bubbles for events/holidays"""
        if special_event == 'christmas':
            return self._create_christmas_bubbles()
        elif special_event == 'halloween':
            return self._create_halloween_bubbles()
        elif special_event == 'valentine':
            return self._create_valentine_bubbles()
        
        return self._create_daily_special_bubbles()
    
    def _create_christmas_bubbles(self):
        """Christmas-themed bubbles"""
        return [
            {'type': 'gift', 'color': '#C62828', 'points': 50, 'effect': 'extra_life'},
            {'type': 'ornament', 'color': '#2E7D32', 'points': 30, 'effect': 'slow_time'},
            {'type': 'candy_cane', 'color': '#FFFFFF', 'points': 40, 'effect': 'double_points'}
        ]
    
    def generate_power_ups(self, game_state):
        """Dynamically generate power-ups based on player needs"""
        power_ups = []
        
        if game_state['lives'] < 2:
            power_ups.append({'type': 'extra_life', 'chance': 0.3})
        
        if game_state['speed_level'] > 8:
            power_ups.append({'type': 'slow_motion', 'chance': 0.4})
        
        if game_state['score'] % 1000 < 100:  # Near milestone
            power_ups.append({'type': 'score_boost', 'chance': 0.5})
        
        return power_ups
